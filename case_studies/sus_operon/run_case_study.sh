#!/usr/bin/env bash
#
# Pynteny case study: the SusC-SusD polysaccharide-utilization pair
# -----------------------------------------------------------------
# Reproduces the whole pipeline end to end:
#   1. download six bacterial genomes (3 Bacteroidota, 3 non-Bacteroidota) from NCBI
#   2. translate + positionally label every ORF with `pynteny build`
#   3. (A) search susC alone, (B) the strict (+)-strand susC-susD tandem,
#      (C) the strand-agnostic tandem, then summarise hit counts.
#
# The two HMMs (PGAP: susC=TIGR04056.1, susD=NF033071.0) and their gene-symbol
# metadata are committed under data/, so the full PGAP database is NOT needed.
#
# Requirements: pynteny (>=1.3.0), curl, gzip, awk.
# Usage:        bash run_case_study.sh
#
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

GENOME_DIR="data/genomes"
PEPTIDE_DIR="data/peptides"
HMM_DIR="data/hmms"
HMM_META="data/sus_hmm_meta.tsv"
RESULTS="results"
NCBI="https://ftp.ncbi.nlm.nih.gov/genomes/all"
# susC/RagA is a large, promiscuous family; the original Sus example uses a permissive
# reporting threshold and lets the synteny filter provide the specificity.
HMMSEARCH_ARGS="-E 1e-10"

mkdir -p "$GENOME_DIR" "$PEPTIDE_DIR" "$RESULTS"

GENOMES=$(awk -F'\t' 'NR>1 && $1!~/^#/ {print $1"\t"$6}' genomes.tsv)

# --- 1 + 2. Download genomes and build labelled peptide databases ----------
while IFS=$'\t' read -r name path; do
    [ -z "$name" ] && continue
    base="$(basename "$path")"
    fna="$GENOME_DIR/${name}.fna"; faa="$PEPTIDE_DIR/${name}.faa"
    if [ ! -s "$fna" ]; then
        echo ">> downloading $name"
        curl -fsSL "$NCBI/$path/${base}_genomic.fna.gz" -o "${fna}.gz"
        gunzip -f "${fna}.gz"
    fi
    if [ ! -s "$faa" ]; then
        echo ">> building labelled peptide database for $name"
        pynteny build --data "$fna" --outfile "$faa"
    fi
done <<< "$GENOMES"

# --- 3. Searches -----------------------------------------------------------
run_search () {           # run_search <label> <outsubdir> <extra-args> <struct>
    local label="$1" sub="$2" extra="$3" struct="$4"
    echo "== $label : '$struct' =="
    while IFS=$'\t' read -r name path; do
        [ -z "$name" ] && continue
        local out="$RESULTS/$sub/$name"
        rm -rf "$out"; mkdir -p "$out"
        # `|| true`: genomes with no hit exit non-zero; that is a valid negative here.
        pynteny search \
            --synteny_struc "$struct" \
            --data "$PEPTIDE_DIR/${name}.faa" \
            --hmm_dir "$HMM_DIR" --hmm_meta "$HMM_META" --gene_ids \
            --hmmsearch_args "$HMMSEARCH_ARGS" \
            --outdir "$out" $extra > "$out/pynteny.log" 2>&1 || true
        local n=0
        [ -s "$out/synteny_matched.tsv" ] && n=$(($(wc -l < "$out/synteny_matched.tsv") - 1))
        printf "   %-40s %s hits\n" "$name" "$n"
    done <<< "$GENOMES"
}

# A) susC alone — a promiscuous TonB-dependent receptor family (SusC/RagA/OmpW),
#    found across many phyla. No genomic context => ambiguous annotation.
run_search "A  susC alone" "A_susC_alone" "" ">susC"

# B) Strict (+)-strand susC-susD tandem — the canonical PUL core, immediately adjacent.
#    Specific to Bacteroidota; catches the ~half of PULs encoded on the (+) strand.
run_search "B  susC-susD tandem, (+)-strand" "B_tandem_strict" "" ">susC 0 >susD"

# C) Strand-agnostic, order-agnostic tandem — recovers the FULL PUL repertoire
#    (both strands), i.e. how many SusC-SusD loci each genome actually carries.
run_search "C  susC-susD tandem, any strand" "C_tandem_any_strand" "--unordered" "susC 0 susD"

# --- 4. Summary ------------------------------------------------------------
SUMMARY="$RESULTS/summary.tsv"
count () { local f="$1"; [ -s "$f" ] && echo $(($(wc -l < "$f") - 1)) || echo 0; }
{
    printf "genome\trole\tsusC_alone\ttandem_strict_members\ttandem_anystrand_members\n"
    while IFS=$'\t' read -r name role; do
        [ -z "$name" ] || [ "$name" = "#name" ] && continue
        printf "%s\t%s\t%s\t%s\t%s\n" "$name" "$role" \
            "$(count "$RESULTS/A_susC_alone/$name/synteny_matched.tsv")" \
            "$(count "$RESULTS/B_tandem_strict/$name/synteny_matched.tsv")" \
            "$(count "$RESULTS/C_tandem_any_strand/$name/synteny_matched.tsv")"
    done < <(awk -F'\t' 'NR>1 {print $1"\t"$2}' genomes.tsv)
} > "$SUMMARY"

echo
echo "Done. Summary written to $SUMMARY:"
column -t -s$'\t' "$SUMMARY"

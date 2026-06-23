#!/usr/bin/env bash
#
# Pynteny case study: detecting the nitrogen-fixation (nif) operon
# ----------------------------------------------------------------
# Reproduces the whole pipeline end to end:
#   1. download six bacterial genomes (3 diazotrophs, 3 non-fixers) from NCBI
#   2. translate + positionally label every ORF with `pynteny build`
#   3. screen each genome for the nitrogenase operon with `pynteny search`
#   4. write per-genome hit tables and a summary into results/
#
# The six nif HMMs (PGAP/TIGRFAM) and their gene-symbol metadata are committed
# under data/, so this script does NOT need the full 432 MB PGAP database.
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
HMM_META="data/nif_hmm_meta.tsv"
RESULTS="results"
NCBI="https://ftp.ncbi.nlm.nih.gov/genomes/all"

mkdir -p "$GENOME_DIR" "$PEPTIDE_DIR" "$RESULTS"

# name <tab> refseq_path  (skip the header and comment lines of genomes.tsv)
GENOMES=$(awk -F'\t' 'NR>1 && $1!~/^#/ {print $1"\t"$6}' genomes.tsv)

# ---------------------------------------------------------------------------
# 1 + 2. Download genomes and build labelled peptide databases
# ---------------------------------------------------------------------------
while IFS=$'\t' read -r name path; do
    [ -z "$name" ] && continue
    base="$(basename "$path")"
    fna="$GENOME_DIR/${name}.fna"
    faa="$PEPTIDE_DIR/${name}.faa"

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

# ---------------------------------------------------------------------------
# 3. Synteny searches
# ---------------------------------------------------------------------------
# --best_hmm_wins is REQUIRED here: nifD/nifK are paralogous to nifE/nifN, so a
# single peptide is often hit by several HMMs; this keeps each peptide's best.

run_search () {           # run_search <label> <outsubdir> <extra-args> <struct>
    local label="$1" sub="$2" extra="$3" struct="$4"
    echo "== $label : '$struct' =="
    while IFS=$'\t' read -r name path; do
        [ -z "$name" ] && continue
        local out="$RESULTS/$sub/$name"
        rm -rf "$out"; mkdir -p "$out"
        # `|| true`: a genome with no matching HMM hits exits non-zero; that is a
        # valid (negative) result for this study, so don't let `set -e` kill us.
        pynteny search \
            --synteny_struc "$struct" \
            --data "$PEPTIDE_DIR/${name}.faa" \
            --hmm_dir "$HMM_DIR" --hmm_meta "$HMM_META" --gene_ids \
            --outdir "$out" $extra > "$out/pynteny.log" 2>&1 || true
        local n=0
        [ -s "$out/synteny_matched.tsv" ] && n=$(($(wc -l < "$out/synteny_matched.tsv") - 1))
        printf "   %-38s %s hits\n" "$name" "$n"
    done <<< "$GENOMES"
}

# A) Strict, strand-aware structural operon: nifH-nifD-nifK, all on the (+)
#    strand, immediately adjacent (0 genes between). Demands exact collinearity.
run_search "A  strict (+)-strand nifHDK" "A_strict_nifHDK" "--best_hmm_wins" \
    ">nifH 0 >nifD 0 >nifK"

# B) Robust, order-/strand-agnostic full panel (structural + cofactor genes).
#    --unordered relaxes both gene order and strand, so it also catches operons
#    on the (-) strand or split by insertion elements (see Nostoc).
run_search "B  unordered nif panel" "B_panel_unordered" "--unordered --best_hmm_wins" \
    "nifB 80 nifH 80 nifD 80 nifK 80 nifE 80 nifN"

# C) Nostoc-specific: the same operon recovered with the correct strand (<) and a
#    wide nifD-nifK gap that spans the ~11-kb nifD excision element.
echo "== C  Nostoc strand/gap-tuned : '<nifK 15 <nifD 0 <nifH' =="
out="$RESULTS/C_nostoc_strand_tuned/Nostoc_PCC7120"
rm -rf "$out"; mkdir -p "$out"
pynteny search --synteny_struc "<nifK 15 <nifD 0 <nifH" \
    --data "$PEPTIDE_DIR/Nostoc_PCC7120.faa" \
    --hmm_dir "$HMM_DIR" --hmm_meta "$HMM_META" --gene_ids \
    --outdir "$out" --best_hmm_wins > "$out/pynteny.log" 2>&1 || true
n=0; [ -s "$out/synteny_matched.tsv" ] && n=$(($(wc -l < "$out/synteny_matched.tsv") - 1))
printf "   %-38s %s hits\n" "Nostoc_PCC7120" "$n"

# ---------------------------------------------------------------------------
# 4. Summary table: genome x search -> number of syntenic hits
# ---------------------------------------------------------------------------
SUMMARY="$RESULTS/summary.tsv"
{
    printf "genome\trole\tA_strict_nifHDK\tB_panel_unordered\n"
    while IFS=$'\t' read -r name role rest; do
        [ -z "$name" ] || [ "$name" = "#name" ] && continue
        a=0; b=0
        f="$RESULTS/A_strict_nifHDK/$name/synteny_matched.tsv"
        [ -s "$f" ] && a=$(($(wc -l < "$f") - 1))
        f="$RESULTS/B_panel_unordered/$name/synteny_matched.tsv"
        [ -s "$f" ] && b=$(($(wc -l < "$f") - 1))
        printf "%s\t%s\t%s\t%s\n" "$name" "$role" "$a" "$b"
    done < <(awk -F'\t' 'NR>1 {print $1"\t"$2}' genomes.tsv)
} > "$SUMMARY"

echo
echo "Done. Summary written to $SUMMARY:"
column -t -s$'\t' "$SUMMARY"

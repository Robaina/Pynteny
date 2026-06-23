# Case study: finding the nitrogen-fixation (*nif*) operon

> *Motivation:* [issue #92](https://github.com/Robaina/Pynteny/issues) — *"This could be great
> for identifying operons. Do you have any examples of using this tool to identify common
> operons? For example, nif operons or others involved in nitrogen fixation."*

This case study uses Pynteny to detect the **nitrogenase (*nif*) operon** — the textbook
example of a syntenic, co-regulated gene cluster — directly from genome assemblies, and to
tell **diazotrophs** (nitrogen fixers) apart from non-fixers.

**Follow along interactively in the notebook: [`nif_operon.ipynb`](nif_operon.ipynb)** — it
walks through every step with the Python API and renders the result tables inline. The six
*nif* HMMs and their metadata are committed under [`data/`](data/), so you do **not** need the
full 432 MB PGAP database.

```bash
conda activate pynteny-dev                 # an env with pynteny >= 1.3.0 + jupyter
jupyter lab nif_operon.ipynb               # interactive walk-through, or
bash run_case_study.sh                      # the same pipeline as a one-shot CLI script
```

---

## 1. The biology

Biological nitrogen fixation — the reduction of atmospheric N₂ to ammonia — is catalysed by
**nitrogenase**, encoded by a compact, strongly conserved gene cluster:

| Gene | HMM (PGAP/TIGRFAM) | Product | Role |
|------|--------------------|---------|------|
| `nifH` | TIGR01287.1 | nitrogenase iron (Fe) protein | electron donor / dinitrogenase reductase |
| `nifD` | TIGR01282.1 | MoFe protein α-chain | catalytic subunit |
| `nifK` | TIGR01286.1 | MoFe protein β-chain | catalytic subunit |
| `nifE` | TIGR01283.1 | NifE | FeMo-cofactor biosynthesis scaffold |
| `nifN` | TIGR01285.1 | NifN | FeMo-cofactor biosynthesis scaffold |
| `nifB` | TIGR01290.1 | NifB | FeMo-cofactor biosynthesis |

The structural genes are almost universally arranged as the **`nifH`–`nifD`–`nifK`** operon,
transcribed as one unit. This co-location and co-orientation is exactly the signal Pynteny is
built to find — far more specific than searching for `nifH`, `nifD` and `nifK` independently.

> ⚠️ **A built-in gotcha — paralogy.** `nifE`/`nifN` are ancient paralogues of `nifD`/`nifK`
> (they arose by duplication of the structural genes). Their HMMs therefore **cross-hit**: one
> peptide is often matched by several models. This study always passes **`--best_hmm_wins`**,
> which keeps each peptide's single highest-scoring HMM and prevents the cross-hits from
> silently breaking the synteny filter.

## 2. The genome panel

Three diazotrophs from three different phyla, and three non-fixers — including one deliberately
counter-intuitive negative (see [`genomes.tsv`](genomes.tsv) for accessions):

| Genome | Role | Lineage | Why it's here |
|--------|------|---------|---------------|
| *Azotobacter vinelandii* DJ | 🟢 fixer | γ-proteobacteria | model free-living aerobic diazotroph; clean `nifHDK` |
| *Sinorhizobium meliloti* 1021 | 🟢 fixer | α-proteobacteria | legume symbiont; *nif* on the **pSymA megaplasmid**, not the chromosome |
| *Nostoc* sp. PCC 7120 | 🟢 fixer | Cyanobacteria | heterocyst-former; `nifD` split by an **11-kb excision element** |
| *Escherichia coli* K-12 MG1655 | 🔴 non-fixer | γ-proteobacteria | classic negative control |
| *Bacillus subtilis* 168 | 🔴 non-fixer | Bacillota | distant-lineage negative control |
| *Klebsiella pneumoniae* 342 | 🔴 non-fixer | γ-proteobacteria | *Klebsiella* is famous for N₂ fixation — **but this strain has no *nif*** |

## 3. Results

Number of syntenic hits per genome ([`results/summary.tsv`](results/summary.tsv)):

| Genome | Role | **A** strict `>nifH 0 >nifD 0 >nifK` | **B** unordered *nif* panel |
|--------|:----:|:-----------------------------------:|:--------------------------:|
| *Azotobacter vinelandii* DJ | 🟢 | **3** | **12** |
| *Sinorhizobium meliloti* 1021 | 🟢 | **3** | **6** |
| *Nostoc* sp. PCC 7120 | 🟢 | 0 | **7** |
| *Escherichia coli* MG1655 | 🔴 | 0 | 0 |
| *Bacillus subtilis* 168 | 🔴 | 0 | 0 |
| *Klebsiella pneumoniae* 342 | 🔴 | 0 | 0 |

**Every diazotroph is detected, every non-fixer is rejected.** The two searches teach
complementary lessons.

### Search A — the strict structural operon

```bash
pynteny search \
    --synteny_struc ">nifH 0 >nifD 0 >nifK" \
    --data data/peptides/<genome>.faa \
    --hmm_dir data/hmms --hmm_meta data/nif_hmm_meta.tsv --gene_ids \
    --best_hmm_wins --outdir results/A_strict_nifHDK/<genome>
```

`>` requires the (+) strand; `0` requires the genes to be immediately adjacent. In
*A. vinelandii* the three genes fall out perfectly collinear
([full table](results/A_strict_nifHDK/Azotobacter_vinelandii_DJ/synteny_matched.tsv)):

| gene № | strand | gene | locus |
|-------:|:------:|------|-------|
| 129 | + | nifH | (136 759, 137 631) |
| 130 | + | nifD | (137 758, 139 236) |
| 131 | + | nifK | (139 337, 140 908) |

*S. meliloti* matches identically — but on contig `NC_003037.1`, the **pSymA megaplasmid**:
Pynteny does not care whether the operon lives on the chromosome or a replicon.

This strict search is **precise but literal**: *Nostoc* scores 0, because its operon is on the
(−) strand *and* is interrupted (next section). That is not a false negative — it is the search
doing exactly what it was told. Which leads to search B.

### Search B — the robust, order- & strand-agnostic panel

```bash
pynteny search \
    --synteny_struc "nifB 80 nifH 80 nifD 80 nifK 80 nifE 80 nifN" \
    --data data/peptides/<genome>.faa \
    --hmm_dir data/hmms --hmm_meta data/nif_hmm_meta.tsv --gene_ids \
    --unordered --best_hmm_wins --outdir results/B_panel_unordered/<genome>
```

No strand symbols and `--unordered` mean "these six genes clustered within ~80 ORFs, in any
order or orientation". This is the recommended **screen** for *"does this genome fix
nitrogen?"*: it recovers all three diazotrophs (and finds *A. vinelandii*'s **two** nitrogenase
copies, 12 hits) while still rejecting every non-fixer.

### The *Nostoc* twist — a 55-year-old discovery, visible in the output

*Nostoc*'s neighbourhood
([full table](results/B_panel_unordered/Nostoc_PCC7120/synteny_matched.tsv)) reads, on the
(−) strand:

| gene № | strand | gene | locus |
|-------:|:------:|------|-------|
| 1476 | − | nifH | (1 713 396, 1 714 283) |
| 1475 | − | nifD | (1 711 821, 1 713 263) |
| … | | **← ~11.5 kb gap, 13 ORFs →** | |
| 1461 | − | nifK | (1 698 743, 1 700 281) |
| 1459 | − | nifE | (1 696 389, 1 697 831) |
| 1458 | − | nifN | (1 695 055, 1 696 389) |
| 1457 | − | nifB | (1 694 529, 1 694 942) |

`nifD` and `nifK` are **not adjacent** — they are separated by ~11.5 kb. This is the famous
**11-kb *nifD* excision element**: in vegetative cells the element interrupts the operon, and
it is excised by the recombinase XisA only during heterocyst differentiation, reconstituting a
functional `nifHDK`. The genome assembly is sequenced from vegetative DNA, so the interruption
is right there in the coordinates.

You can still recover it as an *ordered* operon by telling Pynteny the truth — (−) strand, and
a gap wide enough to span the element (search C in the script):

```bash
pynteny search --synteny_struc "<nifK 15 <nifD 0 <nifH" \
    --data data/peptides/Nostoc_PCC7120.faa \
    --hmm_dir data/hmms --hmm_meta data/nif_hmm_meta.tsv --gene_ids --best_hmm_wins ...
# -> nifK (g1461), nifD (g1475), nifH (g1476), all (-) strand: 3 hits
```

## 4. Takeaways

- **Synteny beats independent gene hits.** Requiring `nifH`–`nifD`–`nifK` *co-located* is far
  more specific than three separate HMM searches — it is the difference between "has a
  nitrogenase-like gene" and "has a nitrogenase operon".
- **Tune three knobs to the biology:** strand (`>` `<` / none), max gene spacing (the integers),
  and order (`--unordered`). Start permissive to screen, tighten to characterise.
- **Mind paralogues.** With cross-hitting models like `nifD`/`nifE` and `nifK`/`nifN`, use
  `--best_hmm_wins`.
- **Presence is strain-specific, not genus-level.** *K. pneumoniae* 342 carries no *nif* even
  though the genus is a textbook nitrogen fixer — only a genome-level search settles it.

## 5. Files

```
nif_operon/
├── README.md                 # this file
├── nif_operon.ipynb          # interactive walk-through (Python API)
├── run_case_study.sh         # the same pipeline as a one-shot CLI script
├── genomes.tsv               # genome panel + NCBI accessions
├── data/
│   ├── hmms/                 # 6 nif HMMs (PGAP/TIGRFAM) — committed (~1.2 MB)
│   ├── nif_hmm_meta.tsv      # HMM accession -> gene symbol map (for --gene_ids)
│   ├── genomes/              # downloaded by the script (git-ignored)
│   └── peptides/             # built by the script   (git-ignored)
└── results/
    ├── summary.tsv           # genome x search hit-count matrix
    ├── A_strict_nifHDK/      # per-genome synteny_matched.tsv + hit FASTAs
    ├── B_panel_unordered/
    └── C_nostoc_strand_tuned/
```

## References & data provenance

- HMMs: NCBI **PGAP** HMM collection (TIGRFAM models), `ftp.ncbi.nlm.nih.gov/hmm/current/`.
- Genomes: NCBI **RefSeq** (accessions in `genomes.tsv`).
- *Nostoc* PCC 7120 *nifD* element: Golden, Robinson & Haselkorn (1985) *Nature* **314**, 419–423.

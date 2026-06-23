# Case study: the SusC-SusD polysaccharide-utilization pair

This case study uses Pynteny to detect **polysaccharide-utilization loci (PULs)**, the gene
clusters that let **Bacteroidota** feed on complex carbohydrates, and to show why **genomic
context beats sequence similarity** for annotating their signature gene, `susC`.

**Follow along interactively in the notebook: [`sus_operon.ipynb`](sus_operon.ipynb)**, it walks
through every step with the Python API and renders the result tables inline. The two HMMs and
their metadata are committed under [`data/`](data/), so you do **not** need the full 432 MB PGAP
database.

```bash
conda activate pynteny-dev                 # an env with pynteny >= 1.3.0 + jupyter
jupyter lab sus_operon.ipynb               # interactive walk-through, or
bash run_case_study.sh                      # the same pipeline as a one-shot CLI script
```

---

## 1. The biology

The pair **`susC`-`susD`** is characteristic of the phylum **Bacteroidota** and lets these
bacteria import long polysaccharides through a "pedal-bin" mechanism. `susC` is a TonB-dependent
outer-membrane transporter; `susD` is the partner substrate-binding lipoprotein. Together they
form the core of a **polysaccharide-utilization locus (PUL)**, and Bacteroidota, abundant in
the ocean, soil and the human gut, often carry **dozens** of them.

The catch: `susC` belongs to a large, widespread family of TonB-dependent receptors (it has
well-known homologs such as **RagA** and **OmpW**), so labelling an unannotated protein "`susC`"
from sequence similarity alone is unreliable. The fix is to use **genomic context**: a `susC`
sitting immediately next to a `susD` is a confident PUL, not just a generic receptor.

| Gene | HMM (PGAP) | Product |
|------|-----------|---------|
| `susC` | TIGR04056.1 | SusC/RagA family TonB-linked outer-membrane protein |
| `susD` | NF033071.0  | starch-binding outer-membrane lipoprotein SusD |

## 2. The genome panel

Three Bacteroidota and three non-Bacteroidota controls (accessions in [`genomes.tsv`](genomes.tsv)):

| Genome | Role | Lineage | Why it's here |
|--------|------|---------|---------------|
| *Bacteroides thetaiotaomicron* VPI-5482 | 🟢 PUL+ | Bacteroidota | the model *Sus* organism; a huge PUL arsenal |
| *Bacteroides fragilis* NCTC 9343 | 🟢 PUL+ | Bacteroidota | another gut symbiont, many PULs |
| *Flavobacterium johnsoniae* UW101 | 🟢 PUL+ | Bacteroidota | environmental (soil/freshwater) Bacteroidota |
| *Pseudomonas aeruginosa* PAO1 | 🔴 control | Pseudomonadota | **many `susC`-like TonB receptors, but no `susD`** |
| *Escherichia coli* K-12 MG1655 | 🔴 control | Pseudomonadota | a few TonB receptors, no `susD` tandem |
| *Bacillus subtilis* 168 | 🔴 control | Bacillota | clean negative (no SusC/SusD system) |

## 3. Results

Hits per genome ([`results/summary.tsv`](results/summary.tsv)):

| Genome | Role | **A** `susC` alone | **B** strict tandem `>susC 0 >susD` (members) | **C** PULs, any strand |
|--------|:----:|:----:|:----:|:----:|
| *B. thetaiotaomicron* VPI-5482 | 🟢 | 69 | 48 | **49** |
| *B. fragilis* NCTC 9343 | 🟢 | 46 | 34 | **37** |
| *F. johnsoniae* UW101 | 🟢 | 32 | 24 | **25** |
| *P. aeruginosa* PAO1 | 🔴 | **6** | 0 | 0 |
| *E. coli* MG1655 | 🔴 | 0 | 0 | 0 |
| *B. subtilis* 168 | 🔴 | 0 | 0 | 0 |

Three searches, three lessons.

### Search A: `susC` alone is promiscuous

```bash
pynteny search --synteny_struc ">susC" --hmmsearch_args "-E 1e-10" \
    --data data/peptides/<genome>.faa \
    --hmm_dir data/hmms --hmm_meta data/sus_hmm_meta.tsv --gene_ids --outdir ...
```

`susC` scores dozens of hits in the Bacteroidota, but also **6 hits in *Pseudomonas
aeruginosa***, whose genome is full of SusC-like TonB receptors that have nothing to do with
PULs. On a lone-`susC` hit you cannot tell a real PUL transporter from a generic receptor.

### Search B: the `susC`-`susD` tandem is specific

```bash
pynteny search --synteny_struc ">susC 0 >susD" --hmmsearch_args "-E 1e-10" ... --outdir ...
```

Requiring `susD` immediately downstream of `susC` (same strand, zero genes between) collapses
the *Pseudomonas* signal to **zero** while keeping the Bacteroidota, turning an ambiguous
family hit into a confident PUL call. In *B. thetaiotaomicron* the hits come out as clean
alternating pairs (`susC` at gene 32 ends at 26 169 bp; `susD` at gene 33 starts at 26 182 bp,
just 13 bp apart).

### Search C: count the whole PUL repertoire

```bash
pynteny search --synteny_struc "susC 0 susD" --unordered --hmmsearch_args "-E 1e-10" ... --outdir ...
```

The strict `(+)`-strand search only sees the ~half of PULs on the forward strand. Dropping the
strand and order constraints recovers them all: **~49 PULs in *B. thetaiotaomicron***, ~37 in
*B. fragilis*, ~25 in *F. johnsoniae*. The gut *Bacteroides*, which must digest a huge diversity
of dietary and host glycans, carry the largest PUL arsenals, exactly as their lifestyle
predicts. Every non-Bacteroidota control stays at zero.

## 4. Takeaways

- **Single-gene hits can mislead.** A lone `susC`-family match (e.g. the six in *Pseudomonas*)
  does not mean "polysaccharide-utilization locus".
- **Synteny supplies the missing specificity.** The adjacent `susC`-`susD` tandem is present in
  the Bacteroidota and **nowhere** in the controls.
- **Tune strand to the question.** `>susC 0 >susD` characterises (+)-strand loci precisely;
  `susC 0 susD --unordered` inventories the entire PUL repertoire, which scales with lifestyle.

## 5. Files

```
sus_operon/
├── README.md                 # this file
├── sus_operon.ipynb          # interactive walk-through (Python API)
├── run_case_study.sh         # the same pipeline as a one-shot CLI script
├── genomes.tsv               # genome panel + NCBI accessions
├── data/
│   ├── hmms/                 # 2 HMMs (susC = TIGR04056.1, susD = NF033071.0), committed
│   ├── sus_hmm_meta.tsv      # HMM accession -> gene symbol map (for --gene_ids)
│   ├── genomes/              # downloaded by the script/notebook (git-ignored)
│   └── peptides/             # built by the script/notebook   (git-ignored)
└── results/
    ├── summary.tsv           # genome x search hit-count matrix
    ├── A_susC_alone/
    ├── B_tandem_strict/
    └── C_tandem_any_strand/
```

## References & data provenance

- HMMs: NCBI **PGAP** collection, `ftp.ncbi.nlm.nih.gov/hmm/current/` (susC = TIGR04056.1,
  susD = NF033071.0).
- Genomes: NCBI **RefSeq** (accessions in `genomes.tsv`).
- This case study generalises the marine-metagenome *Sus* example from the
  [Pynteny docs](https://robaina.github.io/Pynteny/examples/example_sus/) to a small, curated,
  fully reproducible genome panel.

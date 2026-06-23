# Pynteny case studies

Self-contained, reproducible worked examples that apply Pynteny to real biological questions on
real genomes. Each case study is a directory with a `README.md`, an **interactive notebook** to
follow along, a `run_case_study.sh` script that reproduces every result from scratch, the
(small) committed inputs it needs, and the expected outputs under `results/`.

These are heavier and more narrative than the quick-start [`examples/`](../examples) notebooks:
the notebooks show *how to call* Pynteny, the case studies show *how to answer a question* with
it: choosing the HMMs, designing the synteny structure, and interpreting the hits.

| Case study | Question | Highlights |
|------------|----------|------------|
| [`nif_operon/`](nif_operon) | Which genomes fix nitrogen, and where is the nitrogenase operon? | strict vs. unordered searches · paralogue handling (`--best_hmm_wins`) · operon on a megaplasmid · a real operon split by an 11-kb excision element |
| [`sus_operon/`](sus_operon) | Which genomes carry polysaccharide-utilization loci, and how many? | why a syntenic pair beats a single-gene hit · a promiscuous gene family disambiguated by context · counting a whole PUL repertoire across strands |

## Running one

```bash
conda activate pynteny-dev          # any env with pynteny >= 1.3.0 + jupyter
cd nif_operon
jupyter lab nif_operon.ipynb        # interactive, or:
bash run_case_study.sh              # one-shot CLI reproduction
```

## Contributing a case study

Good case studies start from a concrete biological question and a genome panel with a known
expected answer (positives **and** negatives). Keep committed inputs small: commit the HMMs
and metadata, and let the notebook/script download genomes and build peptide databases into
git-ignored `data/` subfolders. See [`nif_operon/`](nif_operon) as the template.

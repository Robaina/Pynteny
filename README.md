![logo](https://user-images.githubusercontent.com/21340147/192824830-dcbe8d09-2b10-431d-bd9a-b4624192dcc9.png)
<br>

# Synteny-aware hmm searches made easy

[![tests](https://github.com/Robaina/Pynteny/actions/workflows/tests.yml/badge.svg)](https://github.com/Robaina/Pynteny/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/Robaina/Pynteny/branch/main/graph/badge.svg?token=WDSOC220X6)](https://codecov.io/gh/Robaina/Pynteny)
[![docs](https://github.com/Robaina/Pynteny/actions/workflows/docs.yml/badge.svg)](https://github.com/Robaina/Pynteny/actions/workflows/docs.yml)

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pynteny/badges/latest_release_date.svg)](https://anaconda.org/bioconda/pynteny)
![license](https://img.shields.io/github/license/Robaina/Pynteny)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)

[![Bioconda](https://img.shields.io/conda/vn/bioconda/pynteny?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pynteny)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pynteny/badges/downloads.svg)](https://anaconda.org/bioconda/pynteny)
[![GitHub release](https://img.shields.io/github/release/Robaina/Pynteny.svg)](https://GitHub.com/Robaina/Pynteny/releases/)


[![Anaconda-Server Badge](https://anaconda.org/bioconda/pynteny/badges/platforms.svg)](https://anaconda.org/bioconda/pynteny)
![python](https://img.shields.io/badge/Python-3.10-blue)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

[![pyOpenSci](https://tinyurl.com/y22nb8up)](https://github.com/pyOpenSci/software-review/issues/67)
[![DOI](https://zenodo.org/badge/500470783.svg)](https://zenodo.org/badge/latestdoi/500470783)


## 1. :bulb: What is Pynteny?

`Pynteny` is Python tool to search for [synteny](https://en.wikipedia.org/wiki/Synteny) blocks in (prokaryotic) sequence data through [HMMs](https://www.bioinformatics.org/wiki/Hidden_Markov_Model) of the ORFs of interest and [HMMER](http://hmmer.janelia.org/). By leveraging genomic context information, `Pynteny` can be employed to decrease the uncertainty of functional annotation of unlabelled sequence data due to the effect of paralogs. `Pynteny` can be accessed (i) through the command line or (ii) as a Python module.

Get more info in the [documentation](https://robaina.github.io/Pynteny/) pages!

## 2. :wrench: Setup

Install with conda:

1. Pynteny requires Python 3.10. The easiest way to handle dependencies is by creating a dedicated conda environment:

```bash
conda create -n pynteny -c bioconda -c conda-forge python=3.10 pynteny
conda activate pynteny
```

2. Check that installation worked fine:

```bash
(pynteny) pynteny --help
```
### 2.1. Installing on Windows

Pynteny is designed to run on Linux machines. However, it can be installed within the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) via conda.

### 2.2. Installing on MacOS with the latest ARM64 architecture

Pynteny doesn't currently support the latest ARM64 architecture of silicon processors (e.g. MacBook M1 and M2). If that is your case, you can install Pynteny using the workaround below (based on [this post](https://towardsdatascience.com/how-to-manage-conda-environments-on-an-apple-silicon-m1-mac-1e29cb3bad12)):

```bash
CONDA_SUBDIR=osx-64 conda create -n pynteny_x86 python=3.10
conda activate pynteny_x86
conda config --env --set subdir osx-64
conda install -c bioconda pynteny
```

## 3. :rocket: Usage

Consider the following toy example of a syntenic block:

![synteny example](assets/synteny_example.png)

Here, we are interested in four genes which colocate according to the pattern above: genes A-C show consecutive locations in the positive strand, followed by three (untargeted) genes and followed by gene D, which is located in the negative strand.

Pynteny can be run either as a command line tool or as a Python module. To run pynteny in the command line, execute:

```bash
conda activate pynteny
pynteny <subcommand> <options>
```

<p align="center">
   <img src="assets/pynteny_cli.png" alt="pynyeny-cli">
</p>


There are a number of available subcommands, which can be explored in the [documentation](https://robaina.github.io/Pynteny/) pages.

For intance, to first download the [PGAP](https://academic.oup.com/nar/article/49/D1/D1020/6018440)'s database containing a collection of profile HMMs as well as metadata:

```bash
pynteny download --outdir data/hmms --unpack
```

Next, to build a labelled peptide database from DNA assembly data:

```bash
pynteny build \
    --data assembly.fa \
    --outfile labelled_peptides.faa

```

Finally, to search the peptide database for the syntenic structure displayed above: `>gene_A 0 >gene_B 0 >gene_C 3 <gene_D`, and using the downloaded [PGAP](https://academic.oup.com/nar/article/49/D1/D1020/6018440) database:

```bash
pynteny search \
    --synteny_struc ">gene_A 0 >gene_B 0 >gene_C 3 <gene_D" \
    --data labelled_peptides.faa \
    --outdir results/ \
    --gene_ids
```

## 4. :notebook_with_decorative_cover: Examples

Here are some Jupyter Notebooks with examples to show how Pynteny works:

<!-- <a href="https://colab.research.google.com/github/Robaina/Pynteny/blob/main/docs/examples/example_api_colab.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a> -->
* [Pynteny API](https://robaina.github.io/Pynteny/examples/example_api/)
* [Pynteny CLI](https://robaina.github.io/Pynteny/examples/example_cli/)
* [Sus operon](https://robaina.github.io/Pynteny/examples/example_sus/)

You can find more notebooks in the [examples directory](docs/examples/). Find more info in the [documentation](https://robaina.github.io/Pynteny/).

## 5. :arrows_counterclockwise: Dependencies
Pynteny would not work without these awesome projects:

- [hmmer](https://github.com/EddyRivasLab/hmmer)
- [prodigal](https://github.com/hyattpd/Prodigal)
- [pyfastx](https://github.com/lmdu/pyfastx)
- [biopython](https://github.com/biopython/biopython)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [numpy](https://github.com/numpy/numpy)
- [pandas](https://github.com/pandas-dev/pandas)
- [psutil](https://github.com/giampaolo/psutil)
- [requests](https://requests.readthedocs.io/en/latest/)
- [tqdm](https://github.com/tqdm/tqdm)

Thanks!

## 6. :octocat: Contributing

Contributions are always welcome! If you don't know where to start, you may find an interesting [issue to work in here](https://github.com/Robaina/Pynteny/issues). Please, read our [contribution guidelines](CONTRIBUTING.md) first.

## 7. :black_nib: Citation

If you use this software, please cite it as below:

Semidán Robaina Estévez. (2023). Pynteny: synteny-aware hmm searches made easy (Version 1.0.0). Zenodo. https://zenodo.org/record/7696204


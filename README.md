![logo](https://user-images.githubusercontent.com/21340147/192824830-dcbe8d09-2b10-431d-bd9a-b4624192dcc9.png)
<br>

# Synteny-aware hmm searches made easy

![python](https://img.shields.io/badge/Python-3.8%2B-blue)
[![tests](https://github.com/Robaina/Pynteny/actions/workflows/tests.yml/badge.svg)](https://github.com/Robaina/Pynteny/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/Robaina/Pynteny/branch/main/graph/badge.svg?token=WDSOC220X6)](https://codecov.io/gh/Robaina/Pynteny)
[![docs](https://github.com/Robaina/Pynteny/actions/workflows/docs.yml/badge.svg)](https://github.com/Robaina/Pynteny/actions/workflows/docs.yml)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pynteny?logo=anaconda&style=flat-square&maxAge=3600)](https://anaconda.org/bioconda/pynteny)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pynteny/badges/downloads.svg)](https://anaconda.org/bioconda/pynteny)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/Robaina/pynteny)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pynteny/badges/platforms.svg)](https://anaconda.org/bioconda/pynteny)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/pynteny/badges/latest_release_date.svg)](https://anaconda.org/bioconda/pynteny)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
![license](https://img.shields.io/github/license/Robaina/Pynteny)
![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4)
[![DOI](https://zenodo.org/badge/500470783.svg)](https://zenodo.org/badge/latestdoi/500470783)
<a href="https://colab.research.google.com/github/Robaina/Pynteny/blob/main/docs/examples/example_api_colab.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>


## 1. :bulb: What is Pynteny?

`Pynteny` is Python tool to search for [synteny](https://en.wikipedia.org/wiki/Synteny) blocks in (prokaryotic) sequence data through [HMMs](https://www.bioinformatics.org/wiki/Hidden_Markov_Model) of the ORFs of interest and [HMMER](http://hmmer.janelia.org/). By leveraging genomic context information, `Pynteny` can be employed to decrease the uncertainty of functional annotation of unlabelled sequence data due to the effect of paralogs. `Pynteny` can be accessed (i) through the command line, (ii) as a Python module or (iii) as a (locally served) web application.

Get more info in the [documentation](https://robaina.github.io/Pynteny/) pages!

## 2. :wrench: Setup

Install with conda:

1. Create dedicated conda environment (highly recommended)

```bash
conda create -n pynteny
conda activate pynteny
```
2. Install pynteny

```bash
conda install -c bioconda pynteny
```

3. Check that installation worked fine:

```bash
pynteny --help
```
### 2.1. Installing Pynteny on Windows

Pynteny is designed to run on Linux machines. However, it can be installed within the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) via conda.

## 3. :rocket: Usage

Pynteny can be run either as a command line tool or as a (locally-served) web application. To run pynteny in the command line, execute:

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

Finally, to search the peptide database for the syntenic structure: `>gene_A 0 >gene_B 0 >gene_C 10 <gene_C`, using the donwloaded [PGAP](https://academic.oup.com/nar/article/49/D1/D1020/6018440) database:

```bash
pynteny search \
    --synteny_struc ">gene_A 0 >gene_B 0 >gene_C 10 <gene_C" \
    --data labelled_peptides.faa \
    --outdir results/ \
    --gene_ids
```

Pynteny may also be used with a graphical interface (made with [Streamlit](https://streamlit.io)). The app is run on a local server in your machine, thus all files are kept locally and the app can be run without an internet connection. 

To run the app, execute the following command once pynteny has been successfully installed:

```bash
conda activate pynteny
pynteny app
```

<p align="center">
   <img src="assets/pynteny_app.gif" alt="pynyeny-app">
</p>

## 4. :notebook_with_decorative_cover: Examples

In the [examples directory](docs/examples/), you can find a collection of Jupyter Notebooks containing workflows to demonstrate the usate of Pynteny's command line interface as well as the Python API. There is also a [Google Colab environment](https://colab.research.google.com/github/Robaina/Pynteny/blob/main/docs/examples/example_api_colab.ipynb) to test Pynteny. We invite you to explore Pynteny's web application by executing the command `pynteny app`. Find more info in the [documentation](https://robaina.github.io/Pynteny/).

## 5. :arrows_counterclockwise: Dependencies
Pynteny would not work without these awesome projects:

- [hmmer](https://github.com/EddyRivasLab/hmmer)
- [prodigal](https://github.com/hyattpd/Prodigal)
- [pyfastx](https://github.com/lmdu/pyfastx)
- [biopython](https://github.com/biopython/biopython)
- [numpy](https://github.com/numpy/numpy)
- [pandas](https://github.com/pandas-dev/pandas)
- [psutil](https://github.com/giampaolo/psutil)
- [python-wget](https://anaconda.org/conda-forge/python-wget/)
- [streamlit](https://github.com/streamlit/streamlit)
- [streamlit-aggrid](https://github.com/PablocFonseca/streamlit-aggrid)

Thanks!

## 6. :octocat: Contributing

Contributions are always welcome! If you don't know where to start, you may find an interesting [issue to work in here](https://github.com/Robaina/Pynteny/issues). Please, read our [contribution guidelines](CONTRIBUTING.md) first.

## 7. :black_nib: Citation

If you use this software, please cite it as below:

Semidán Robaina Estévez. (2022). Pynteny: synteny-aware hmm searches made easy (Version 0.0.4). Zenodo. https://doi.org/10.5281/zenodo.7048685


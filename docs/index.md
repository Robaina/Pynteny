![logo](https://user-images.githubusercontent.com/21340147/191948930-3ee11b4e-4b13-4365-b14a-2b709c76a49f.png)

Pynteny is a Python package to query sequence databases for synteny blocks through Hidden Markov Models (HMMs). It can be run in a command line interface or through a graphical interface. The command line tool is more complete and is intended to be employed in remote servers. In contrast, the graphical interface runs locally and may be easier to handle for some users (also intended for educational purposes).

These are the available subcommands, run as ```pynteny <subcommand> <options>```:

- [search](subcommands/search.md)
- [build](subcommands/build.md)
- [download](subcommands/download.md)
- [parse](subcommands/parse.md)
- [app](subcommands/app.md)

## Setup

Install with conda:

1. Create dedicated conda environment (highly recommended)

```bash
conda create -n pynteny
conda activate pynteny
```
2. Install pynteny

Pynteny can be install with the following command:
```bash
conda install -c bioconda pynteny
```
If conda takes a long time to solve the dependencies (e.g., more than a few minutes) or fails to solve it, you may try:
```bash
conda install -c conda-forge -c bioconda python pynteny
```
Or install it with [mamba](https://github.com/mamba-org/mamba) instead of conda:
```bash
mamba install -c bioconda pynteny
```

3. Check that installation worked fine:

```bash
pynteny --help
```
### Installing Pynteny on Windows

Pynteny is designed to run on Linux machines. However, it can be installed within the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install) via conda.

## General usage

Pynteny's main subcommand, [`pynteny search`](subcommands/search.md), requires a peptide (ORF) sequence database in fasta format in which record labels contain positional information of each sequence with respect to their contig of origin. Additionally, labels must follow the following format:

```
<genome ID>_<contig ID>_<gene position>_<locus start>_<locus end>_<strand>
```

To facilitate the construction of this database, Pynteny provides the subcommand [`pynteny build`](subcommands/build.md), which takes as input a fasta file containing assembled nucleotide sequence data (or a single or a collection of genomes), such as that retrieved from MAG reconstruction pipelines.

Additionally, [`pynteny search`](subcommands/search.md) requires either the set of profile Hidden Markov Models (HMMs) used in the provided synteny structure or a database of profile HMMs from which to retrieve the necessary HMMs. The user may provide their own set of HMMs. However, you may also download the entire [PGAP HMM database](https://academic.oup.com/nar/article/49/D1/D1020/6018440) through [`pynteny download`](subcommands/download.md), which will take care of downloading and storing the download path for future usage (by default any time no additional HMMs are provided as an argument in [`pynteny search`](subcommands/search.md)).

Once both the peptide database and the required HMMs are ready, you can query the peptide database with a text string encoding the query synteny structure such as the following:

$>HMM_a \:\: n_{ab} \:\: > (HMM_{b1} | HMM_{b2}) \:\: n_{bc} \:\: < HMM_c.$

Where $HMM_a$ represents the name of the HMM to be used (corresponding to the file name without the extension), $n_{ab}$ is an integer representing the maximum number of genes between HMMs a and b, < and > indicate the strand in which to search for the HMM pattern, antisense and sense, respectively. Note that more than one HMM can be employed for a single gene in the structure, as indicated by the HMM group $(HMM_{b1} | HMM_{b2})$ above. In these cases, [`pynteny search`](subcommands/search.md) will search for sequences that matched any HMM contained within the HMM group.

## Getting started and Examples

Here are some Jupyter Notebooks with examples to show how Pynteny works.

<a href="https://colab.research.google.com/github/Robaina/Pynteny/blob/main/docs/examples/example_api_colab.ipynb" target="_blank"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>
* [Pynteny API](https://robaina.github.io/Pynteny/examples/example_api/)
* [Pynteny CLI](https://robaina.github.io/Pynteny/examples/example_cli/)
* [Sus operon](https://robaina.github.io/Pynteny/examples/example_sus/)

You can find more notebooks in the [examples directory](docs/examples/). We invite you to explore Pynteny's web application by executing the command `pynteny app`. Find more info in the [documentation](https://robaina.github.io/Pynteny/).

## Contributing

Contributions are always welcome! If you don't know where to start, you may find an interesting [issue to work in here](https://github.com/Robaina/Pynteny/issues). Please, read our [contribution guidelines](https://github.com/Robaina/Pynteny/blob/main/CONTRIBUTING.md) first.

## Citation

If you use this software, please cite it as below:

Semidán Robaina Estévez. (2022). Pynteny: synteny-aware hmm searches made easy(Version 0.0.4). Zenodo. https://doi.org/10.5281/zenodo.7048685
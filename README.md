<img style="width: 100vw; height: 250px" src="assets/pynteny_logo.png">
<br>

# Synteny hmm searches made easy
Pynteny is a command-line tool to search for [synteny](https://en.wikipedia.org/wiki/Synteny) blocks in (assembled) sequence data through [HMMs](https://www.bioinformatics.org/wiki/Hidden_Markov_Model) of the ORFs of interest and [HMMER](http://hmmer.janelia.org/).

## Installation

Download (or fork) repo to local directory:

```bash
git clone https://github.com/Robaina/Pynteny.git
```

cd to downloaded repo and install conda enviroment:

```bash
cd Pynteny
conda env create -f environment.yml
```

Install pynteny in conda enviroment:

```bash
conda activate pynteny
python setup.py install
```

Check that installation worked fine:

```bash
pynteny --help
```

## Usage

Pynteny currently contains two main subcommands:

* `pynteny search`: searches for synteny blocks in a set of ORFs using HMMER and outputs the results in a tabular format. Synteny blocks are specified by strings of ordered HMM names or gene IDs with the following format:

```
$>HMM_a n_{ab} >HMM_b n_{bc} <HMM_c$
```

* `pynteny preprocess`: predict ORFs with [prodigal]()and add positional information to each ORF &mdash; i.e., loci and gene number within assembled contig. Alternatively, the user can provide their own ORF annotations in GenBank format.
* `pynteny download`: downloads the latest version of the [NCBI Prokaryotic Genome Annotation Pipeline](https://github.com/ncbi/pgap) (PGAP) HMM database. However, the user may provide their own HMM database.



<br>
<br>
<br>
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

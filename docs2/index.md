![logo](https://user-images.githubusercontent.com/21340147/191948930-3ee11b4e-4b13-4365-b14a-2b709c76a49f.png)

Pynteny is a Python package to query sequence databases for synteny blocks through Hidden Markov Models (HMMs). It can be run in a command line interface or through a graphical interface. The command line tool is more complete and is intended to be employed in remote servers. In contrast, the graphical interface runs locally and may be easier to handle for some users (also intended for educational purposes).

These are the available subcommands, run as ```pynteny <subcommand> <options>```:

- [search](https://github.com/Robaina/Pynteny/wiki/search)
- [build](https://github.com/Robaina/Pynteny/wiki/build)
- [download](https://github.com/Robaina/Pynteny/wiki/download)
- [parse](https://github.com/Robaina/Pynteny/wiki/parse)
- [tests](https://github.com/Robaina/Pynteny/wiki/tests)
- [app](https://github.com/Robaina/Pynteny/wiki/app)

## General usage

Pynteny's main subcommand, [`pynteny search`](https://github.com/Robaina/Pynteny/wiki/search), requires a peptide (ORF) sequence database in fasta format in which record labels contain positional information of each sequence with respect to their contig of origin. Additionally, labels must follow the following format:

```
<genome ID>_<contig ID>_<gene position>_<locus start>_<locus end>_<strand>
```

To facilitate the construction of this database, Pynteny provides the subcommand [`pynteny build`](https://github.com/Robaina/Pynteny/wiki/build), which takes as input a fasta file containing assembled nucleotide sequence data (or a single or a collection of genomes), such as that retrieved from MAG reconstruction pipelines.

Additionally, [`pynteny search`](https://github.com/Robaina/Pynteny/wiki/search) requires either the set of profile Hidden Markov Models (HMMs) used in the provided synteny structure or a database of profile HMMs from which to retrieve the necessary HMMs. The user may provide their own set of HMMs. However, the user may also download the entire [PGAP HMM database](https://academic.oup.com/nar/article/49/D1/D1020/6018440) through [`pynteny download`](https://github.com/Robaina/Pynteny/wiki/download), which will take care of downloading and storing the download path for future usage (by default any time no additional HMMs are provided as an argument in [`pynteny search`](https://github.com/Robaina/Pynteny/wiki/search)).

Once both the peptide database and the required HMMs are ready, the user can query the peptide database with a text string encoding the query synteny structure such as the following:

$>HMM_a \\:\\: n_{ab} \\:\\: > (HMM_{b1} | HMM_{b2}) \\:\\: n_{bc} \\:\\: < HMM_c.$

Where $HMM_a$ represents the name of the HMM to be used (corresponding to the file name without the extension), $n_{ab}$ is an integer representing the maximum number of genes between HMMs a and b, < and > indicate the strand in which to search for the HMM pattern, antisense and sense, respectively. Note that more than one HMM can be employed for a single gene in the structure, as indicated by the HMM group $(HMM_{b1} | HMM_{b2})$ above. In these cases, [`pynteny search`](https://github.com/Robaina/Pynteny/wiki/search) will search for sequences that matched any HMM contained within the HMM group.

## Examples

There are some Jupyter notebooks with examples of how to use Pynteny's command-line interface as well as its Python API. Check them out [here](https://github.com/Robaina/Pynteny/tree/master/examples)!
## Citation

If you use this software, please cite it as below:

Semidán Robaina Estévez. (2022). Pynteny: synteny-aware hmm searches made easy(Version 0.0.2). Zenodo. https://doi.org/10.5281/zenodo.7048685
Translate assembly sequence data and assign positional metadata to labels

## Usage

```bash
usage: pynteny build [-h] [args] 

```
## Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
|`-i`|`--data`|`None`|path to assembly input nucleotide data or annotated GenBank file.  It can be a single file or a directory of files (either of FASTA or GeneBank format).|
|`-o`|`--outfile`|`None`|path to output (labeled peptide database) file. Defaults to file in the directory of input data.|
|`-n`|`--processes`|`None`|set the maximum number of processes. Defaults to all but one.|
|`-l`|`--log`|`None`|path to log file. Log not written by default.|


## Description

Translate nucleotide assembly file and assign contig and gene location info 
to each identified ORF (using [prodigal](https://github.com/hyattpd/Prodigal)). Label predicted ORFs according to 
positional info and export a fasta file containing predicted and translated ORFs. 
Alternatively, extract peptide sequences from GenBank file containing ORF annotations 
and write labeled peptide sequences to a fasta file.
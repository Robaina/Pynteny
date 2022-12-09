Query sequence database for HMM hits arranged in provided synteny structure.

## Usage

```bash
usage: pynteny search [-h] [args] 

```
## Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
|`-s`|`--synteny_struc`|`None`|string displaying hmm structure to search for, such as:    '>hmm_a n_ab <hmm_b n_bc hmm_c'   where '>' indicates a hmm target located on the positive strand,  '<' a target located on the negative strand, and n_ab corresponds to the maximum number of genes separating matched genes a and b.  Multiple hmms may be employed.  No order symbol in a hmm indicates that results should be independent of strand location. |
|`-i`|`--data`|`None`|path to peptide database|
|`-d`|`--hmm_dir`|`None`|path to directory containing hmm (i.e, tigrfam or pfam) models.  The directory can contain more hmm models than used in the synteny structure.  It may also be the path to a compressed (tar, tar.gz, tgz) directory.  If not provided, hmm models (PGAP database) will be downloaded from the NCBI. (if not already downloaded)|
|`-o`|`--outdir`|`None`|path to output directory|
|`-x`|`--prefix`|``|prefix to be added to output files|
|`-p`|`--processes`|`None`|maximum number of processes available to HMMER. Defaults to all but one.|
|`-a`|`--hmmsearch_args`|`None`|list of comma-separated additional arguments to hmmsearch for each input hmm.  A single argument may be provided, in which case the same additional argument is employed in all hmms.|
|`-g`|`--gene_ids`||use gene symbols in synteny structure instead of HMM names.  If set, a path to the hmm database metadata file must be provided in argument '--hmm_meta'|
|`-u`|`--unordered`||whether the HMMs should be arranged in the exact same order displayed in the synteny_structure or in any order. If ordered, the filters will filter collinear rather than syntenic structures. If more than two HMMs are employed, the largest maximum distance among any pair is considered to run the search.|
|`-r`|`--reuse`||reuse hmmsearch result table in following synteny searches.  Do not delete hmmer_outputs subdirectory for this option to work.|
|`-m`|`--hmm_meta`|`None`|path to hmm database metadata file|
|`-l`|`--log`|`None`|path to log file. Log not written by default.|


## Description

Search for synteny blocks in a set of ORFs using HMMER and outputs the results in a tabular format. Synteny blocks are specified by strings of ordered HMM names or gene IDs with the following format:

$$\lt HMM_a \space n_{ab} \space \lt HMM_b \space n_{bc} \space \lt(HMM_{c1}|HMM_{c2}|HMM_{c3}),$$ 

where $n_{ab}$ corresponds to the maximum number of genes between $HMM_a$ and $HMM_b$. Results can be strand-specific, in that case, $>$ preceding an HMM name indicates that the corresponding ORF must be located in the positive (or sense) strand. Likewise, a $<$ symbol indicates that the ORF must be located in the negative (antisense) strand. Searches can be made strand-insensitive by omitting the $>$ or $<$ symbol. 

Several HMMs can be assigned to the same ORF, in which case the search is performed for all of them. In this case, HMM names must be separated by "|" and grouped within parentheses, as shown above.

If the PGAP database is employed (see `pynteny download` below), synteny blocks can also be specified by gene symbols, such as $$\lt leuD \space 0 \space \lt leuC \space 1 \space \lt leuA.$$ In that case, the program will try to match gene symbols to HMM names in the PGAP database before running the search.

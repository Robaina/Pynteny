Download PGAP's HMM database from NCBI.

# Usage:

```bash
usage: pynteny download [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
|`-o`|`--outdir`|`None`|path to the directory where to download HMM database.|
|`-u`|`--unpack`||unpack originally compressed database files|
|`-f`|`--force` ||force-download database again if already downloaded.|
|`-l`|`--log`|`None`|path to log file. Log not written by default.|

## Description

download the latest version of the [NCBI Prokaryotic Genome Annotation Pipeline](https://github.com/ncbi/pgap) (PGAP) HMM database. Users may provide their own HMM database in the [search](https://github.com/Robaina/Pynteny/wiki/search) subcommand.
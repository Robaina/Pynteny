# Description

Translate synteny structure with gene symbols into one with
HMM groups, according to provided HMM database.

# Usage:

```bash
usage: pynteny parse [-h] [args] 

```
# Arguments

|short|long|default|help|
| :--- | :--- | :--- | :--- |
|`-h`|`--help`||show this help message and exit|
|`-s`|`--synteny_struc`|`None`|synteny structure containing gene symbols instead of HMMs|
|`-m`|`--hmm_meta`|`None`|path to hmm database metadata file. If already downloaded with pynteny download, hmm meta file is retrieved from the default location.|
|`-l`|`--log`|`None`|path to log file. Log not written by default.|
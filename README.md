# Pynteny: synteny hmm searches made easy in Python

1. Install environment

conda env create -f environment.yml

# Possible applications:

1. Detect operons in (assembled) environmental samples, MAGs for instance. These operons may be specific of certain taxa, which could lead to taxonomic identification/confirmation besides functional identification.

# Current limitations:
1. If using several, alternative HMMs (for the same gene), pynteny only retains matches


# TODO: 
2. Enable search of TigrFAM by gene name? need gene/TigrFAM dictionary
3. Use python logger instead of print statements

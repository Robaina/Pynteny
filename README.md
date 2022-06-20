# Pynteny: synteny hmm searches made easy in Python

1. Install environment

conda env create -f environment.yml

# Possible applications:

1. Detect operons in (assembled) environmental samples, MAGs for instance. These operons may be specific of certain taxa, which could lead to taxonomic identification/confirmation besides functional identification.


# TODO: 
1. Enable search by several hmms for the same gene, like: >(A|B) 0 < C:
   1.1 Need to modify removing duplicated hits, since two hmms for the same gene will produce duplicated hits
   1.2 Probably best way is to preprocess hits from equivalent hmms. Then treat these hits as a single "combined" gene
2. Enable search of TigrFAM by gene name? need gene/TigrFAM dictionary
3. Organize code into a package / CLI
4. Use python logger instead of print statements

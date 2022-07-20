# Example with E. coli MG1655
# https://biocyc.org/ECOLI/NEW-IMAGE?type=LOCUS-POSITION&object=NIL&chromosome=COLI-K12&orgids=ECOLI&bp-range=70001/120000
#
# leuA -> TIGR00973.1
# leuB -> TIGR00169.1
# leuC -> TIGR00170.1
# leuD -> TIGR00171.1
# setA -> TIGR00899.1
# dapF -> TIGR00652.1
#
# >setA 0 <leuD 0 <leuC 0 <leuB 0 <leuA

# rm tests/MG1655/*

pynteny download --unpack

pynteny search \
 --data pynteny/tests/data/MG1655.fasta \
 --outdir pynteny/tests/MG1655/ \
 --synteny_struc "<leuD 0 <leuC 1 <leuA" --gene_ids
 
 # The command should retrieve the following synteny block:
 # "<(TIGR00171.1|TIGR02084.1) 0 <(TIGR00170.1|TIGR02083.1) 1 <(NF002084.0|TIGR00973.1|TIGR00970.1)"
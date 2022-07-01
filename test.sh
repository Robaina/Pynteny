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

# rm tests/MG1655_results/*

# pynteny search \
#  --hmm_dir /home/robaina/Databases/hmm_PGAP/ \
#  --in /home/robaina/Documents/Pynteny/tests/MG1655.fasta \
#  --outdir /home/robaina/Documents/Pynteny/tests/MG1655_results \
#  --synteny_struc "<leuD 0 <leuC 1 <leuA" --gene_ids --hmm_meta hmm_PGAP.tsv     #  --synteny_struc "<(TIGR00171.1|TIGR02084.1) 0 <(TIGR00170.1|TIGR02083.1) 1 <(NF002084.0|TIGR00973.1|TIGR00970.1)"
 

# Searching MAR database (SAR11 synteny structure)
# <leuB 0 <leuD 0 <leuC 4 >dapF

# python translate_assembly.py \
#  --assembly_fasta /home/robaina/Databases/MAR_database/marref_assembly_V6_example.fa \
#  --outdir /home/robaina/Databases/MAR_database/ \
#  --prefix "marref_V6_" \
#  --processes 7 \
#  --split_contigs --metagenome

# python pynteny.pynteny.py preprocess \
#  --assembly_fasta /home/robaina/Databases/MAR_database/marref_assembly_V6_example.fa \
#  --outdir /home/robaina/Databases/MAR_database/ \
#  --prefix "marref_V6_" \
#  --processes 7 \
#  --split_contigs --metagenome

# python synteny_search.py \
#  --hmm_dir /home/robaina/Databases/hmm_PGAP/ \
#  --in /home/robaina/Databases/MAR_database/marref_V6_positioned_clean.faa \
#  --outdir /home/robaina/Documents/Pynteny/tests/sar11_results \
#  --prefix "sar11_" \
#  --synteny_struc "<TIGR00170.1 10 >TIGR00652.1"  # "<TIGR00169.1 0 <TIGR00171.1 0 <TIGR00170.1"

# Searching for sox operon in MAR database 
# >soxX 0 >soxY 0 >soxZ 0 >soxA

# rm tests/sox_results/*

# pynteny search \
#  --in /home/robaina/Databases/MAR_database/marref_prodigal_longlabels.faa \
#  --outdir /home/robaina/Documents/Pynteny/tests/sox_results_2 \
#  --synteny_struc ">soxX 0 >soxY 0 >soxZ 0 >soxA 0 >soxB 0 >soxC" \
#  --gene_ids  # --hmm_meta /home/robaina/Databases/hmm_PGAP.tsv \
#  #  --hmm_dir /home/robaina/Databases/hmm_PGAP/ \




pynteny search \
 --in /home/robaina/Databases/MAR_database/marref_prodigal_longlabels.faa \
 --outdir /home/robaina/Documents/Pynteny/tests/sox_results_2 \
 --synteny_struc ">soxX 0 >soxY 0 >soxZ 0 >soxA 0 >soxB 0 >soxC" \
 --gene_ids
#  --hmm_dir /home/robaina/Documents/Pynteny/data/hmm_PGAP.HMM.tgz \
#  --hmm_meta /home/robaina/Documents/Pynteny/data/hmm_PGAP.tsv \
#  --gene_ids
 
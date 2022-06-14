# Pynteny: synteny hmm searches made easy in Python

1. Genome annotation provided via gbk/gff3 file, via text file with contig/gene info or via prokka if desired and prokaryotic organism

# Possible applications:

1. Detect operons in (assembled) environmental samples, MAGs for instance. These operons may be specific of certain taxa, which could lead to taxonomic identification/confirmation besides functional identification.


# TODO: 
1. Redesign filtering strategy
2. Use Path from pathlib instead of str

# Packages:
1. pandas
2. pyfastx
3. biopython
4. prodigal
5. hmmer
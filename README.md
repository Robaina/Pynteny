# Pynteny: synteny hmm searches made easy in Python

1. Install environment

conda env create -f environment.yml

# Possible applications:

0. Detect operons in (assembled) environmental samples, MAGs for instance. These operons may be specific of certain taxa, which could lead to taxonomic identification/confirmation besides functional identification.

A roadmap for metagenomic enzyme discovery

https://pubs.rsc.org/en/content/articlelanding/2021/np/d1np00006c


1. Biosynthetic Gene Clusters (BGCs):Minimum Information about a BGC (MIBiG) database (https://academic.oup.com/nar/article/48/D1/D454/5587631)

2. sequence similarity networks (SSNs) are relatively new methods for the visualization of protein families and superfamilies. First published for the purpose of protein superfamily analysis in 2009,131 SSNs are graphs that display relationships between protein families. 

3. to automate genome neighborhood analysis, a widely used addition to the EFI-EST is the Genome Neighborhood Tool (GNT)

4. Based on its ‘BiG’ savings in computational cost, BiG-SLICE is therefore particularly well-suited for analysis of metagenomes for genome-context guided enzyme discovery.

5. Taking this one step further, deep learning methods for embedding genes as vectors in their genomic context (e.g., pfam2vec) have led to improvements in BGC prediction.


# Current limitations:
1. If using several, alternative HMMs (for the same gene), pynteny only retains matches


# TODO: 
2. Double-check final output, some of the hmms in each hmm group may not have returned any hit yet they appear in results (HMM1 | HMM2)
3. Use python logger instead of print statements

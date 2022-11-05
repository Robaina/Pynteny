---
title: 'Pynteny: a Python package to perform synteny-aware, profile HMM-based searches in sequence databases'
tags:
  - Python
  - bioinformatics
  - HMMER
  - synteny
  - HMMs
  - sequencing
authors:
  - name: Semidán Robaina-Estévez
    orcid: 0000-0003-0781-1677
    email: srobaina@ull.edu.es
    affiliation: 1
affiliations:
 - name: Department of Microbiology. University of La Laguna. Spain.
   index: 1
date: 3 November 2022
bibliography: paper.bib
---


# Summary

With an exponentially growing number of available sequence data, automated function annotation of sequences has become a key subfield of Bioinformatics. In a majority of cases, annotation methods rely on sequence similarity to sequences with known functions to assign functional labels &mdash; assuming that similarity implies shared ancestry, that is, homology. Sequence similarity is most commonly assessed by either alignment-based methods, such as BLAST [@blast], or sequence profile-based methods, such as HMMER3 [@hmmer]. In the first case, query sequences are aligned and compared to a reference sequence database. In the latter, however, query sequences are compared to a profile Hidden Markov Model (HMM), a probabilistic model of the sequence space which is obtained from a collection of representative sequences with the same annotated function. Therefore, profile-based methods are particularly well-suited when query sequences are not well represented in reference databases, as it facilitates the search of distant homologs due to the natural variability encoded in the profile HMM [@hmms-a,@hmms-b].

While function is generally conserved among sequence orthologs, i.e., homologs that are the result of a speciation event, this is not the general case of paralogs, that is, homologs that are the result of a gene duplication event, which typically undergo functional diversification. Due to the existence of paralogs, it is impossible to assess orthology solely based on sequence similarity, and additional sources of information, such as phylogenetics and genomic context are necessary to resolve paralogous from orthologous sequences. The consideration of genomic context, such as _synteny_ &mdash; the physical co-location of genes within the same chromosome across different species &mdash; during function annotation is particularly useful in prokaryotes, where genes tend to cluster together into operons. In these cases, syntenic information can reduce annotation uncertainty by providing additional, co-localization constraints to the homology search.

Constraining profile-based searches with syntenic information could markedly benefit annotation pipelines of prokaryotic sequences, particularly those originating in metagenomic samples which typically are poorly represented in reference databases.

# Statement of need

`Pynteny` is a Python tool to conduct synteny-aware, profile HMM searches in prokaryotic sequence databases. To this end, it allows querying sequence databases with arrangements of profile HMMs that reflect a target syntenic block instead of with single profile HMMs. `Pynteny` is designed to work directly with assembled sequence data, it relies on Prodigal [@prodigal] to translate and add positional tags to individual genes, and on HMMER3 [@hmmer] to search sequence databases for homologs through profile HMMs, providing Python wrappers for both programs.

`Pynteny` was designed to be used by researchers working with large, unannotated sequence databases, such as those typically encountered in metagenomic analyses. Researchers who feel comfortable with Python may integrate `Pynteny` in their pipelines through its API. Those who do not feel comfortable with programming may use `Pynteny` through its command line or web interface, although the latter presents limited functionality and it is more suitable for educational purposes.

# State of the field

Several existing tools are dedicated to the exploration, analysis, and visualization of synteny blocks among genomes. In these tools, users typically input a number of annotated genomes and obtain a collection of syntenic relations of shared gene sets among the genomes. Examples of these tools are MCScan [@mcscan] and MCScanX [@mcscanx], Clinker [@clinker], pyGenomeViz [@pyGenomeViz], genePlotR [@geneplotR], gggenomes [@gggenomes], GENESPACE [@genespace], Mology [@mology]. While these tools are excellent resources to identify and analyze syntenic relations among genomes, they are functionally complementary to `Pynteny`. Specifically, rather than exploring syntenic blocks within annotated genomes, `Pynteny`'s objective is to search for specific syntenic structures within unannotated (assembled) sequence data, as well as to leverage syntenic information to reduce uncertainty due to paralogs during function annotation. Therefore, `Pynteny` requires a previous identification of conserved syntenic structures which can be done with the previously indicated tools.

# Acknowledgements

We acknowledge constructive feedback from Pynteny users that has helped improve the package. This work has been conducted without any financial or commercial support.

# References
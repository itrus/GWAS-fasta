# GWAS
R code for creating Manhattan plots for GWAS

This code is made for analysis of short sequences (0.001-1000 kb) presented in FASTA format. It needs two files (plus.fas and minus.fas) for making SNP analysis of mutations associated with phenotypical difference observed in positive (plus.fas) and control samples (minus.fas). It utilizes Fisher's or Chi squared (default) test.

Use gwas_aa.R for aminoacid sequences and gwas_nt.R for nucleotide ones. Speed of work is approximately 1 Mb of fasta files (plus.fas+minus.fas) in 1 minute on Core 2 Duo 2.66 GHz/4G RAM.

Required packages: seqinr, ade4.

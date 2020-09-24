# Leptidea_hostplant_project
Scripts used in the Leptidea hostplant-project for differential gene expression and microbiome analysis. 
The pipeline for RNA-seq processing was created by Luis Leal (ref) and the microbiome analysis was created by Axel Künstner (ref).

The workflow is described in the manuscript but these are the main steps. 

# RNA-seq processing and differential gene expression analysis

FASTQC version 0.11.5 (Andrews, 2016) - overall quality assesment of the RNA-seq reads.
TRIMGALORE version 0.4.4 (Krueger, 2017) - initial filtering of raw RNA-seq reads including trimming 12 nucleotide bases from the 5’ end of each sequence, trimming sequences with overall Phred score < 30, filtering out adapter sequences and sequences shorter than 30 bp. 
FASTQ MASKER (http://hannonlab.cshl.edu/fastx_toolkit/; accessed 2019-05-01) - mask (replace with N) low quality (threshold = 10) nucleotides in the reads. 
PRINSEQ version 0.20.4 (Schmieder & Edwards, 2011) - trim remaining poly-A tails.
CUTADAPT version 2.5 (Martin, 2011) - filter out long stretches of A/T nucleotides (threshold = 10bp) inside reads. 
CONDETRI (Smeds & Künstner, 2011) - filter out reads with Phred Score < 30 in more than 80% of the read. 
SortMeRNA version 2.1 (Kopylova et al., 2012) - removing remaining sequences of ribosomal origin (rRNA).
FASTQ Screen (Wingett, 2017) with Bowtie2 version 2.3.5 (Langmead & Salzberg, 2012) - identify and screen for contaminants. The most likely contaminants were identified based on previous gene expression studies in L. sinapis (Höök et al., 2019; Leal et al., 2018) including: rRNA, human, Drosophila melanogaster, Wolbachia sp., L. corniculatus, L. dorycnium, L. japonicus, Illumina adapters and primers.
STAR version 2.7.2b (Dobin et al., 2013) - indexing and mapping using previously available L. sinapis genome (Talla et al., 2017) and transcriptome (Höök et al., 2019; Leal et al., 2018) assemblies. 
STRINGTIE version 1.3.6 (Pertea et al., 2015) - obtain gene specific counts of the mapped reads.  
DESeq2 version 3.6 (Love et al., 2014) as implemented in R version 3.4.3 (R Core Team, 2013) - standardized differential gene expression analyses using the raw gene counts. Genes for which only 1 sample in the entire sample set had non-zero read counts and genes with zero counts in all samples in a specific cohort were removed. A count of one (1) was added to every gene/sample, to stabilize variance at low expressed genes. Low coverage genes with baseMean (count average across all samples) < 2 (to account for the 1 count added to every gene/sample) were removed before carrying out the differential expression analysis. DESeq2 was run with default settings. This protocol normalizes counts per gene by library size (the number of reads in a specific library) and carries out significance testing for individual genes using the Wald test (Love et al., 2014). The analysis implements the method of Benjamini and Hochberg (1995) to generates false discovery rate (FDR) adjusted significance levels (padj) for each gene. 
topGO version 2.38.1 (Alexa & Rahnenfuhrer, 2016) in the Bioconductor package in R version 3.4.3 (R Core Team, 2013) - assess enrichment of specific gene functions in differentially expressed gene sets using the database org.Dm.eg.db from the R-package AnnotationDbi (Pagès et al., 2019) for orthologous genes for the gene ontology categories biological process, cellular component and molecular function.

# Microbiome characterisation

dada2 version 1.14.0 (Callahan et al., 2016) -  read filtering and processing. Omitting low quality reads, clipping reads to 295 bp (maxEE = 1, maxN = 0, truncQ = 2 and phiX removal activated) and removing chimeric sequences with removeBimeraDenovo function in dada2. 
IDTaxa (Murali et al., 2018) in the R package DECIPHER version 2.14.0 (Wright, 2016) - taxonomy assignment to each amplicon sequence variant (ASV) using GTDB (r89; Parks et al., 2018) as reference database. 
AlignSeqs from the DECIPHER package and phanghorn version 2.5.5 (Schliep, 2011) - alignment of the ASV:s and phylogenetic tree reconstruction applying a general time reversible model with gamma optimization and stochastic rearrangement (optim.pml command). ASVs with kingdom not belonging to bacteria, or unknown phylum, or belonging to family Mitochondria or genus Chloroplast were removed. 
phyloseq version 1.30.0 (McMurdie & Holmes, 2013) and lme package version 3.1-143) - estimate alpha diversity (Shannon, 1948) and evaluate cohort effects in a linear mixed model framework and non-parametric test (Wilcoxon-test). Beta diversity was assessed using the unweighted UniFrac (function in phyloseq) distance (Lozupone & Knight, 2005) on subsampled data (rarefied to 500 contigs per sample), as suggested by Weiss et al. (2017). Permutational multivariate analysis of variance using distance matrices (PERMANOVA) was performed using the adonis command in vegan version 2.5-6 (Oksanen et al., 2019) with 99,999 bootstrap permutations. Differences in taxonomic abundances were investigated using non-parametric testing (Wilcoxon’s test). All analysis were performed in R version 3.6.1 and p-values were corrected using Benjamini-Hochberg correction (p.adjust function in R).

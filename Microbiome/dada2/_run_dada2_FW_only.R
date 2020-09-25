# https://benjjneb.github.io/dada2/

# install package
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.9")

library(dada2)
library(phyloseq)

library(DECIPHER)
library(phangorn)

#
# Filter data
#

# File parsing
path <- "00_fastq/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

plotQualityProfile(fnFs[1:4])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "../01_filtered_FW", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, truncLen=c(295),
                     maxN=0, maxEE=c(1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE, verbose = TRUE, randomize = TRUE, MAX_CONSIST = 20)

plotErrors(errF, nominalQ=TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE, verbose = TRUE)

# Construct sequence table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 440:466]

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)

WriteXLS::WriteXLS(data.frame(track), ExcelFileName = "dada2.fwonly.track.xlsx", row.names = TRUE, col.names = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1)


# Assign taxonomy
# tax <- assignTaxonomy(seqtab.nochim, "~/Databases/Microbiome/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# using decipher
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
#load("~/Databases/Microbiome/SILVA_SSU_r132_March2018.RData") 
load("~/Databases/Microbiome/GTDB_r89-mod_June2019.RData") 
ids <- IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
    m <- match(ranks, x$rank)
    taxa <- x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
}))
colnames(taxid) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus","Species") 
rownames(taxid) <- getSequences(seqtab.nochim)

head(taxid)

# Construct Phylogenetic Tree
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences
alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit <- pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


# Write to disk
saveRDS(seqtab.nochim, "FW_final_seqtab.rds") # CHANGE ME to where you want sequence table saved
saveRDS(taxid, "FW_final_tax.rds") # CHANGE ME ...
saveRDS(fitGTR, "FW_final_tree.rds") # CHANGE ME ...



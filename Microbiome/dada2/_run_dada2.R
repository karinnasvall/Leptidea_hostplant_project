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
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2)

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:4])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "../01_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "../01_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(290,240),
                     maxN=0, maxEE=c(2,3), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
#head(mergers[[1]])

# Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 440:466]

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab2)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

WriteXLS::WriteXLS(data.frame(track), ExcelFileName = "dada2.track.xlsx", row.names = TRUE, col.names = TRUE, BoldHeaderRow = TRUE, FreezeRow = 1)


# Assign taxonomy
# tax <- assignTaxonomy(seqtab.nochim, "~/Databases/Microbiome/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
# using decipher
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("~/Databases/Microbiome/SILVA_SSU_r132_March2018.RData")
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
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))


# Write to disk
saveRDS(seqtab.nochim, "final_seqtab.rds") # CHANGE ME to where you want sequence table saved
saveRDS(taxid, "final_tax.rds") # CHANGE ME ...
saveRDS(fitGTR, "final_tree.rds") # CHANGE ME ...










############
#
# Additional stuff....
#


ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
             #  sample_data(samdf),
               tax_table(taxid),
             phy_tree(fitGTR$tree))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

plot_richness(ps, measures=c("Shannon", "Chao1"))





# remve lowly abundant ASVs
x <- subset_taxa(ps, Kingdom=="Bacteria")
x <- subset_taxa(x, Phylum!="unclassified")
x <- subset_taxa(x, Phylum!="Bacteria_unclassified")
sort( apply ( data.frame(otu_table(x)), 1, sum ) )
fil <- genefilter_sample( x, filterfun_sample(function(x) x > 1), A = 5 ); sum(fil)
( x <- prune_taxa(fil, x) )

sort( apply ( data.frame(otu_table(x)), 1, sum ) )

x <- rarefy_even_depth(physeq = x, sample.size = 9000, rngseed = 123, trimOTUs = TRUE, replace = FALSE)
x

otus <- t(data.frame(otu_table(x)))
otus <- data.frame( "OTU" = rownames(otus) , otus )
rownames(otus) <- NULL

asv.phylum <- t( data.frame( otu_table( tax_glom(physeq = x, taxrank = "Phylum") ) ) )
asv.family <- t( data.frame( otu_table( tax_glom(physeq = x, taxrank = "Family") ) ) )
asv.genera <- t( data.frame( otu_table( tax_glom(physeq = x, taxrank = "Genus")  ) ) )

rownames( asv.phylum ) <- c( tax_table( tax_glom(physeq = x, taxrank = "Phylum") )[, "Phylum"] )
rownames( asv.family ) <- c( tax_table( tax_glom(physeq = x, taxrank = "Family") )[, "Family"] )
rownames( asv.genera ) <- c( tax_table( tax_glom(physeq = x, taxrank = "Genus")  )[, "Genus"] )

write.table( x = otus, file = "asv_table.csv", quote = F, sep = "\t", row.names = F )
write.table( x = asv.phylum, file = "phylum_table.csv", quote = F, sep = "\t", row.names = T )
write.table( x = asv.family, file = "family_table.csv", quote = F, sep = "\t", row.names = T )
write.table( x = asv.genera, file = "genus_table.csv",  quote = F, sep = "\t", row.names = T )

Biostrings::writeXStringSet(refseq(x), 'asv_seqs.fasta')



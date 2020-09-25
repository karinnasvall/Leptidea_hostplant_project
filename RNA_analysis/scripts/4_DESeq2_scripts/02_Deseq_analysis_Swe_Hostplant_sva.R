# Differential gene expression using DEseq2

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
biocLite("sva")
biocLite("cowplot")
library("DESeq2")
library(sva)


sessionInfo()

#Clear all states
rm(list=ls(all=TRUE))
dev.off()


#Load gene count matrix and labels ####

setwd("/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/4_DESeq2")
countData <- as.matrix(read.csv("/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/2c_Quantification_StringTie_STEP2/gene_count_matrix.csv", row.names="gene_id"))
head(countData)
nrow(countData)

#setwd("/home/luisleal/MYPROJ/3_DosageCompensation_LS/Scripts_NEW")
colData <- read.csv("/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/scripts/4_DESeq2_scripts/sample_info_DEseq.txt", sep="\t", row.names=1, header=F)
colnames(colData) <- c("Hostplant", "Population", "Stage", "Family")
colData

#Note: The PHENO_DATA file contains information on each sample, e.g., sex or population.
#The exact way to import this depends on the format of the file.




##### output folders

# standard filtering
OUT_FOLD_STDF <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/4_DESeq2/05_Swe_Hostplant_effect/"

# minimal filtering
#OUT_FOLD_MINF <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/NEW_09_DESEQ2_sva_minFiltering"



#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####  >>> NEW (in the original script this was being done individually for each stage)
countData <- countData[rowSums(countData > 0) >= 2 , ]
nrow(countData)



#Make separate data sets for each developmental stage ####
Instar3_Swe_Host_countdata <- subset(countData, select = c("P11917_143_S90", "P11917_146_S93", "P11917_149_S96", "P11917_152_S99", "P11917_155_S102", 
                                                  "P11917_158_S105", "P11917_161_S108", "P11917_164_S111"))
Instar3_Swe_Host_coldata <- as.data.frame(colData[c(42, 45, 48, 51, 54, 57, 60, 63), ])

Instar5M_Swe_Host_countdata <- subset(countData, select = c("P11917_144_S91", "P11917_147_S94", "P11917_150_S97", "P11917_153_S100", "P11917_156_S103", 
                                                            "P11917_159_S106", "P11917_162_S109", "P11917_165_S112"))

Instar5M_Swe_Host_coldata <- as.data.frame(colData[c(43, 46, 49, 52, 55, 58, 61, 64), ])

Instar5F_Swe_Host_countdata <- subset(countData, select = c("P11917_145_S92", "P11917_148_S95", "P11917_151_S98", "P11917_154_S101", "P11917_157_S104", 
                                                            "P11917_160_S107", "P11917_163_S110", "P11917_166_S113"))
Instar5F_Swe_Host_coldata <- as.data.frame(colData[c(44, 47, 50, 53, 56, 59, 62, 65), ])

Instar3_Swe_Host_coldata <- Instar3_Swe_Host_coldata[ , 1:4]
Instar5M_Swe_Host_coldata <- Instar5M_Swe_Host_coldata[ , 1:4]
Instar5F_Swe_Host_coldata <- Instar5F_Swe_Host_coldata[ , 1:4]



#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####    >>> this is now done above, for all samples (as opposed to doing it for each the individual stage)
#Instar_countdata <- Instar_countdata[rowSums(Instar_countdata > 0) >= 2 , ]
#Pupa_countdata <- Pupa_countdata[rowSums(Pupa_countdata > 0) >= 2 , ]
#Adult_countdata <- Adult_countdata[rowSums(Adult_countdata > 0) >= 2 , ]

#Pre-filtering II: for each developmental stage, remove genes with zero counts in all samples ####  
Instar3_Swe_Host_countdata <- Instar3_Swe_Host_countdata[rowSums(Instar3_Swe_Host_countdata > 0) >= 1 , ]
nrow(Instar3_Swe_Host_countdata)
Instar5M_Swe_Host_countdata <- Instar5M_Swe_Host_countdata[rowSums(Instar5M_Swe_Host_countdata > 0) >= 1 , ]
nrow(Instar5M_Swe_Host_countdata)
Instar5F_Swe_Host_countdata <- Instar5F_Swe_Host_countdata[rowSums(Instar5F_Swe_Host_countdata > 0) >= 1 , ]
nrow(Instar5F_Swe_Host_countdata)

## Add 1 count to every gene/sample (multiple zero-count samples take DESEQ off track)
Instar3_Swe_Host_countdata <- Instar3_Swe_Host_countdata + 1
head(Instar3_Swe_Host_countdata)
Instar5M_Swe_Host_countdata <- Instar5M_Swe_Host_countdata + 1
head(Instar5M_Swe_Host_countdata)
Instar5F_Swe_Host_countdata <- Instar5F_Swe_Host_countdata + 1
head(Instar5F_Swe_Host_countdata)





#Check all sample IDs in colData are also in CountData and match their orders ####
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

all(rownames(Instar3_Swe_Host_coldata) %in% colnames(Instar3_Swe_Host_countdata))
Instar3_Swe_Host_countdata <- Instar3_Swe_Host_countdata[, rownames(Instar3_Swe_Host_coldata)]
all(rownames(Instar3_Swe_Host_coldata) == colnames(Instar3_Swe_Host_countdata))

all(rownames(Instar5M_Swe_Host_coldata) %in% colnames(Instar5M_Swe_Host_countdata))
Instar5M_Swe_Host_countdata <- Instar5M_Swe_Host_countdata[, rownames(Instar5M_Swe_Host_coldata)]
all(rownames(Instar5M_Swe_Host_coldata) == colnames(Instar5M_Swe_Host_countdata))

all(rownames(Instar5F_Swe_Host_coldata) %in% colnames(Instar5F_Swe_Host_countdata))
Instar5F_Swe_Host_countdata <- Instar5F_Swe_Host_countdata[, rownames(Instar5F_Swe_Host_coldata)]
all(rownames(Instar5F_Swe_Host_coldata) == colnames(Instar5F_Swe_Host_countdata))




#Create a DESeqDataSet from count matrix and labels (no confounding effects removed at this stage) ####
ddsInstar3_Swe_Host <- DESeqDataSetFromMatrix(countData = Instar3_Swe_Host_countdata,
                                                 colData = Instar3_Swe_Host_coldata, design = ~Hostplant)
ddsInstar3_Swe_Host

ddsInstar5M_Swe_Host <- DESeqDataSetFromMatrix(countData = Instar5M_Swe_Host_countdata,
                                                  colData = Instar5M_Swe_Host_coldata, design = ~Hostplant)
ddsInstar5M_Swe_Host

ddsInstar5F_Swe_Host <- DESeqDataSetFromMatrix(countData = Instar5F_Swe_Host_countdata,
                                                  colData = Instar5F_Swe_Host_coldata, design = ~Hostplant)
ddsInstar5F_Swe_Host



#Pre-filtering IIIs: Only keep rows with baseMean >2 (usually this woulf be baseMean >1 but we added one count to every gene/sample)
ddsInstar3_Swe_Host <- ddsInstar3_Swe_Host[rowMeans(counts(ddsInstar3_Swe_Host)) > 2, ]
ddsInstar3_Swe_Host
ddsInstar5M_Swe_Host <- ddsInstar5M_Swe_Host[rowMeans(counts(ddsInstar5M_Swe_Host)) > 2, ]
ddsInstar5M_Swe_Host
ddsInstar5F_Swe_Host <- ddsInstar5F_Swe_Host[rowMeans(counts(ddsInstar5F_Swe_Host)) > 2, ]
ddsInstar5F_Swe_Host




######### Run the default analysis for DESeq2 ####

# DESeq2 (Without batch effect correction)

ddsInstar3_Swe_Host <- DESeq(ddsInstar3_Swe_Host)
res_Instar3_Swe_Host <- results(ddsInstar3_Swe_Host, alpha=.05)
res_Instar3_Swe_Host

ddsInstar5M_Swe_Host <- DESeq(ddsInstar5M_Swe_Host)
res_Instar5M_Swe_Host <- results(ddsInstar5M_Swe_Host, alpha=.05)
res_Instar5M_Swe_Host

ddsInstar5F_Swe_Host <- DESeq(ddsInstar5F_Swe_Host)
res_Instar5F_Swe_Host <- results(ddsInstar5F_Swe_Host, alpha=.05)
res_Instar5F_Swe_Host






### Filtering ###


#Filter baseMean > 10 ####
Filtered_base_Mean_10_Instar3_Swe_Host <- subset(res_Instar3_Swe_Host, baseMean>10)
Filtered_base_Mean_10_Instar5M_Swe_Host <- subset(res_Instar5M_Swe_Host, baseMean>10)
Filtered_base_Mean_10_Instar5F_Swe_Host <- subset(res_Instar5F_Swe_Host, baseMean>10)

summary(Filtered_base_Mean_10_Instar3_Swe_Host)
summary(Filtered_base_Mean_10_Instar5M_Swe_Host)
summary(Filtered_base_Mean_10_Instar5F_Swe_Host)


# Filter logFC > |1.0|) ####
Filtered_base_Mean_10_LFC_1.0_Instar3_Swe_Host <- subset(Filtered_base_Mean_10_Instar3_Swe_Host, log2FoldChange >1 | log2FoldChange < -1)
Filtered_base_Mean_10_LFC_1.0_Instar5M_Swe_Host <- subset(Filtered_base_Mean_10_Instar5M_Swe_Host, log2FoldChange >1 | log2FoldChange < -1)
Filtered_base_Mean_10_LFC_1.0_Instar5F_Swe_Host <- subset(Filtered_base_Mean_10_Instar5F_Swe_Host, log2FoldChange >1 | log2FoldChange < -1)

summary(Filtered_base_Mean_10_LFC_1.0_Instar3_Swe_Host)
summary(Filtered_base_Mean_10_LFC_1.0_Instar5M_Swe_Host)
summary(Filtered_base_Mean_10_LFC_1.0_Instar5F_Swe_Host)

## Filter p-adj < 0.05 (no batch effects correction)
Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Swe_Host <- subset(Filtered_base_Mean_10_LFC_1.0_Instar3_Swe_Host, padj<0.05)
Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Swe_Host <- subset(Filtered_base_Mean_10_LFC_1.0_Instar5M_Swe_Host, padj<0.05)
Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Swe_Host <- subset(Filtered_base_Mean_10_LFC_1.0_Instar5F_Swe_Host, padj<0.05)

summary(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Swe_Host)
summary(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Swe_Host)
summary(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Swe_Host)



#######




###### Write results to file ####


######### standard filtering (genes with p-adj < 0.05, |logFC| > 1, base_Mean > 10)

setwd(OUT_FOLD_STDF)

#unfiltered
write.table(res_Instar3_Swe_Host,
            file = "res_Instar3_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(res_Instar5M_Swe_Host, 
            file = "res_Instar5M_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(res_Instar5F_Swe_Host, 
            file = "res_Instar5F_Swe_Host.txt", sep = "\t", quote = FALSE)


#filtered baseMean |1|
write.table(Filtered_base_Mean_10_Instar3_Swe_Host, 
            file = "Filtered_base_Mean_10_Instar3_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_Instar5M_Swe_Host, 
            file = "Filtered_base_Mean_10_Instar5M_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_Instar5F_Swe_Host, 
            file = "Filtered_base_Mean_10_Instar5F_Swe_Host.txt", sep = "\t", quote = FALSE)


#filtered LFC |1|
write.table(Filtered_base_Mean_10_LFC_1.0_Instar3_Swe_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_Instar3_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_Instar5M_Swe_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_Instar5M_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_Instar5F_Swe_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_Instar5F_Swe_Host.txt", sep = "\t", quote = FALSE)


#filtered p<0.05
write.table(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Swe_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Swe_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Swe_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Swe_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Swe_Host.txt", sep = "\t", quote = FALSE)






##########


######### 




### Plots ###


install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')

library(EnhancedVolcano)


# Volcano Plot

EnhancedVolcano(res_Instar3_SpaSwe_Host,
                lab = rownames(res_Instar3_SpaSwe_Host),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8))


EnhancedVolcano(Filtered_base_Mean_10_Instar5F_Swe_Host,
                lab = rownames(Filtered_base_Mean_10_Instar5F_Swe_Host),
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-8, 8),
                title = 'Filtered_base_Mean_10_Instar5F_Swe_Host',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0)


### Plotting Adjusted p values

EnhancedVolcano(Filtered_base_Mean_10_Instar5F_Swe_Host,
                lab = rownames(Filtered_base_Mean_10_Instar5F_Swe_Host),
                x = 'log2FoldChange',
                y = 'padj',
                xlim=c(-10,10),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                title = 'Swe_Instar5F',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                transcriptLabSize = 0,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)


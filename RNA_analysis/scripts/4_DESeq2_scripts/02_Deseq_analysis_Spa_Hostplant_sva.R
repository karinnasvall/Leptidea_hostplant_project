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
OUT_FOLD_STDF <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/4_DESeq2/06_Spa_Hostplant_effect/"

# minimal filtering
#OUT_FOLD_MINF <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/NEW_09_DESEQ2_sva_minFiltering"



#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####  >>> NEW (in the original script this was being done individually for each stage)
countData <- countData[rowSums(countData > 0) >= 2 , ]
nrow(countData)



#Make separate data sets for each developmental stage ####
Instar3_Spa_Host_countdata <- subset(countData, select = c("P11917_101_S49", "P11917_104_S52", "P11917_107_S55", 
                                                  "P11917_110_S58", "P11917_113_S61", "P11917_116_S64", "P11917_119_S67", "P11917_125_S72", 
                                                  "P11917_128_S75", "P11917_131_S78", "P11917_134_S81", "P11917_137_S84", "P11917_140_S87"))
Instar3_Spa_Host_coldata <- as.data.frame(colData[c(1, 4, 7, 10, 13, 16, 19, 24, 27, 30, 33, 36, 39), ])

Instar5M_Spa_Host_countdata <- subset(countData, select = c("P11917_102_S50", "P11917_105_S53", "P11917_108_S56", "P11917_111_S59", "P11917_114_S62", 
                                                   "P11917_117_S65", "P11917_120_S68", "P11917_123_S70", "P11917_126_S73", "P11917_129_S76", 
                                                   "P11917_132_S79", "P11917_135_S82", "P11917_138_S85", "P11917_141_S88"))

Instar5M_Spa_Host_coldata <- as.data.frame(colData[c(2, 5, 8, 11, 14, 17, 20, 22, 25, 28, 31, 34, 37, 40), ])

Instar5F_Spa_Host_countdata <- subset(countData, select = c("P11917_103_S51", "P11917_106_S54", "P11917_109_S57", "P11917_112_S60", "P11917_115_S63", 
                                                   "P11917_118_S66", "P11917_121_S69", "P11917_124_S71", "P11917_127_S74", "P11917_130_S77", 
                                                   "P11917_133_S80", "P11917_136_S83", "P11917_139_S86", "P11917_142_S89"))
Instar5F_Spa_Host_coldata <- as.data.frame(colData[c(3, 6, 9, 12, 15, 18, 21, 23, 26, 29, 32, 35, 38, 41), ])

Instar3_Spa_Host_coldata <- Instar3_Spa_Host_coldata[ , 1:4]
Instar5M_Spa_Host_coldata <- Instar5M_Spa_Host_coldata[ , 1:4]
Instar5F_Spa_Host_coldata <- Instar5F_Spa_Host_coldata[ , 1:4]



#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####    >>> this is now done above, for all samples (as opposed to doing it for each the individual stage)
#Instar_countdata <- Instar_countdata[rowSums(Instar_countdata > 0) >= 2 , ]
#Pupa_countdata <- Pupa_countdata[rowSums(Pupa_countdata > 0) >= 2 , ]
#Adult_countdata <- Adult_countdata[rowSums(Adult_countdata > 0) >= 2 , ]

#Pre-filtering II: for each developmental stage, remove genes with zero counts in all samples ####  
Instar3_Spa_Host_countdata <- Instar3_Spa_Host_countdata[rowSums(Instar3_Spa_Host_countdata > 0) >= 1 , ]
nrow(Instar3_Spa_Host_countdata)
Instar5M_Spa_Host_countdata <- Instar5M_Spa_Host_countdata[rowSums(Instar5M_Spa_Host_countdata > 0) >= 1 , ]
nrow(Instar5M_Spa_Host_countdata)
Instar5F_Spa_Host_countdata <- Instar5F_Spa_Host_countdata[rowSums(Instar5F_Spa_Host_countdata > 0) >= 1 , ]
nrow(Instar5F_Spa_Host_countdata)

## Add 1 count to every gene/sample (multiple zero-count samples take DESEQ off track)
Instar3_Spa_Host_countdata <- Instar3_Spa_Host_countdata + 1
head(Instar3_Spa_Host_countdata)
Instar5M_Spa_Host_countdata <- Instar5M_Spa_Host_countdata + 1
head(Instar5M_Spa_Host_countdata)
Instar5F_Spa_Host_countdata <- Instar5F_Spa_Host_countdata + 1
head(Instar5F_Spa_Host_countdata)






#Check all sample IDs in colData are also in CountData and match their orders ####
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

all(rownames(Instar3_Spa_Host_coldata) %in% colnames(Instar3_Spa_Host_countdata))
Instar3_Spa_Host_countdata <- Instar3_Spa_Host_countdata[, rownames(Instar3_Spa_Host_coldata)]
all(rownames(Instar3_Spa_Host_coldata) == colnames(Instar3_Spa_Host_countdata))

all(rownames(Instar5M_Spa_Host_coldata) %in% colnames(Instar5M_Spa_Host_countdata))
Instar5M_Spa_Host_countdata <- Instar5M_Spa_Host_countdata[, rownames(Instar5M_Spa_Host_coldata)]
all(rownames(Instar5M_Spa_Host_coldata) == colnames(Instar5M_Spa_Host_countdata))

all(rownames(Instar5F_Spa_Host_coldata) %in% colnames(Instar5F_Spa_Host_countdata))
Instar5F_Spa_Host_countdata <- Instar5F_Spa_Host_countdata[, rownames(Instar5F_Spa_Host_coldata)]
all(rownames(Instar5F_Spa_Host_coldata) == colnames(Instar5F_Spa_Host_countdata))





#Create a DESeqDataSet from count matrix and labels (no confounding effects removed at this stage) ####
ddsInstar3_Spa_Host <- DESeqDataSetFromMatrix(countData = Instar3_Spa_Host_countdata,
                                              colData = Instar3_Spa_Host_coldata, design = ~Hostplant)
ddsInstar3_Spa_Host

ddsInstar5M_Spa_Host <- DESeqDataSetFromMatrix(countData = Instar5M_Spa_Host_countdata,
                                               colData = Instar5M_Spa_Host_coldata, design = ~Hostplant)
ddsInstar5M_Spa_Host

ddsInstar5F_Spa_Host <- DESeqDataSetFromMatrix(countData = Instar5F_Spa_Host_countdata,
                                               colData = Instar5F_Spa_Host_coldata, design = ~Hostplant)
ddsInstar5F_Spa_Host


#Pre-filtering IIIs: Only keep rows with baseMean >2 (usually this woulf be baseMean >1 but we added one count to every gene/sample)
ddsInstar3_Spa_Host <- ddsInstar3_Spa_Host[rowMeans(counts(ddsInstar3_Spa_Host)) > 2, ]
ddsInstar3_Spa_Host
ddsInstar5M_Spa_Host <- ddsInstar5M_Spa_Host[rowMeans(counts(ddsInstar5M_Spa_Host)) > 2, ]
ddsInstar5M_Spa_Host
ddsInstar5F_Spa_Host <- ddsInstar5F_Spa_Host[rowMeans(counts(ddsInstar5F_Spa_Host)) > 2, ]
ddsInstar5F_Spa_Host



######### Run the default analysis for DESeq2 ####

# DESeq2 (Without batch effect correction)

ddsInstar3_Spa_Host <- DESeq(ddsInstar3_Spa_Host)
res_Instar3_Spa_Host <- results(ddsInstar3_Spa_Host, alpha=.05)
res_Instar3_Spa_Host

ddsInstar5M_Spa_Host <- DESeq(ddsInstar5M_Spa_Host)
res_Instar5M_Spa_Host <- results(ddsInstar5M_Spa_Host, alpha=.05)
res_Instar5M_Spa_Host

ddsInstar5F_Spa_Host <- DESeq(ddsInstar5F_Spa_Host)
res_Instar5F_Spa_Host <- results(ddsInstar5F_Spa_Host, alpha=.05)
res_Instar5F_Spa_Host





### Filtering ###


#Filter baseMean > 10 ####
Filtered_base_Mean_10_Instar3_Spa_Host <- subset(res_Instar3_Spa_Host, baseMean>10)
Filtered_base_Mean_10_Instar5M_Spa_Host <- subset(res_Instar5M_Spa_Host, baseMean>10)
Filtered_base_Mean_10_Instar5F_Spa_Host <- subset(res_Instar5F_Spa_Host, baseMean>10)

summary(Filtered_base_Mean_10_Instar3_Spa_Host)
summary(Filtered_base_Mean_10_Instar5M_Spa_Host)
summary(Filtered_base_Mean_10_Instar5F_Spa_Host)


# Filter logFC > |1.0|) ####
Filtered_base_Mean_10_LFC_1.0_Instar3_Spa_Host <- subset(Filtered_base_Mean_10_Instar3_Spa_Host, log2FoldChange >1 | log2FoldChange < -1)
Filtered_base_Mean_10_LFC_1.0_Instar5M_Spa_Host <- subset(Filtered_base_Mean_10_Instar5M_Spa_Host, log2FoldChange >1 | log2FoldChange < -1)
Filtered_base_Mean_10_LFC_1.0_Instar5F_Spa_Host <- subset(Filtered_base_Mean_10_Instar5F_Spa_Host, log2FoldChange >1 | log2FoldChange < -1)

summary(Filtered_base_Mean_10_LFC_1.0_Instar3_Spa_Host)
summary(Filtered_base_Mean_10_LFC_1.0_Instar5M_Spa_Host)
summary(Filtered_base_Mean_10_LFC_1.0_Instar5F_Spa_Host)


# Filter p-adj < 0.05 (no batch effects correction)
Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Spa_Host <- subset(Filtered_base_Mean_10_LFC_1.0_Instar3_Spa_Host, padj<0.05)
Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Spa_Host <- subset(Filtered_base_Mean_10_LFC_1.0_Instar5M_Spa_Host, padj<0.05)
Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Spa_Host <- subset(Filtered_base_Mean_10_LFC_1.0_Instar5F_Spa_Host, padj<0.05)

summary(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Spa_Host)
summary(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Spa_Host)
summary(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Spa_Host)



#######



###### Write results to file ####


######### standard filtering (genes with p-adj < 0.05, |logFC| > 1, base_Mean > 10)

setwd(OUT_FOLD_STDF)

#unfiltered
write.table(res_Instar3_Spa_Host,
            file = "res_Instar3_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(res_Instar5M_Spa_Host, 
            file = "res_Instar5M_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(res_Instar5F_Spa_Host, 
            file = "res_Instar5F_Spa_Host.txt", sep = "\t", quote = FALSE)


#filtered baseMean |1|
write.table(Filtered_base_Mean_10_Instar3_Spa_Host, 
            file = "Filtered_base_Mean_10_Instar3_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_Instar5M_Spa_Host, 
            file = "Filtered_base_Mean_10_Instar5M_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_Instar5F_Spa_Host, 
            file = "Filtered_base_Mean_10_Instar5F_Spa_Host.txt", sep = "\t", quote = FALSE)

#filtered LFC |1|
write.table(Filtered_base_Mean_10_LFC_1.0_Instar3_Spa_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_Instar3_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_Instar5M_Spa_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_Instar5M_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_Instar5F_Spa_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_Instar5F_Spa_Host.txt", sep = "\t", quote = FALSE)


#filtered p<0.05
write.table(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Spa_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Spa_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Spa_Host.txt", sep = "\t", quote = FALSE)
write.table(Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Spa_Host, 
            file = "Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5F_Spa_Host.txt", sep = "\t", quote = FALSE)






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

EnhancedVolcano(Filtered_base_Mean_10_Instar5F_Spa_Host,
                lab = rownames(Filtered_base_Mean_10_Instar5F_Spa_Host),
                x = 'log2FoldChange',
                y = 'padj',
                xlim=c(-10,10),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                title = 'Spa_Instar5F',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                transcriptLabSize = 0,
                legend=c('NS','Log2 FC','Adjusted p-value',
                         'Adjusted p-value & Log2 FC'),
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0)


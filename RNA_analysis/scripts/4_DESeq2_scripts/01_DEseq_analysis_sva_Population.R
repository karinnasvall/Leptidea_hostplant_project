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

setwd("/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/4_DESeq2/02_Population_effect")
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
OUT_FOLD_STDF <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Veronika/RNA_HostPlantProject/results/4_DESeq2/02_Population_effect/"

# minimal filtering
#OUT_FOLD_MINF <- "/crex1/proj/uppstore2017185/b2014034_nobackup/Luis/3_DosageCompensation_LS/NEW_09_DESEQ2_sva_minFiltering"




#From deseq2 manual:
#...for some datasets, exploratory data analysis (EDA) plots could reveal that one
#or more groups has much higher within-group variability than the others. A simulated
#example of such a set of samples is shown below. This is case where, by comparing 
#groups A and B separately - subsetting a DESeqDataSet to only samples from those two 
#groups and then running DESeq on this subset - will be more sensitive than a model 
#including all samples together. It should be noted that such an extreme range of 
#within-group variability is not common, although it could arise if certain treatments
#produce an extreme reaction (e.g. cell death). Again, this can be easily detected
#from the EDA plots such as PCA described in this vignette.





#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####  >>> NEW (in the original script this was being done individually for each stage)
countData <- countData[rowSums(countData > 0) >= 2 , ]
nrow(countData)



#Make separate data sets for each developmental stage ####
Instar3_countdata <- subset(countData, select = c("P11917_101_S49", "P11917_104_S52", "P11917_107_S55", 
                                                 "P11917_110_S58", "P11917_113_S61", "P11917_116_S64", "P11917_119_S67", "P11917_125_S72", 
                                                 "P11917_128_S75", "P11917_131_S78", "P11917_134_S81", "P11917_137_S84", "P11917_140_S87", 
                                                 "P11917_143_S90", "P11917_146_S93", "P11917_149_S96", "P11917_152_S99", "P11917_155_S102", 
                                                 "P11917_158_S105", "P11917_161_S108", "P11917_164_S111"))
Instar3_coldata <- as.data.frame(colData[c(1, 4, 7, 10, 13, 16, 19, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63), ])

Instar5M_countdata <- subset(countData, select = c("P11917_102_S50", "P11917_105_S53", "P11917_108_S56", "P11917_111_S59", "P11917_114_S62", 
                                                   "P11917_117_S65", "P11917_120_S68", "P11917_123_S70", "P11917_126_S73", "P11917_129_S76", 
                                                   "P11917_132_S79", "P11917_135_S82", "P11917_138_S85", "P11917_141_S88", "P11917_144_S91", 
                                                   "P11917_147_S94", "P11917_150_S97", "P11917_153_S100", "P11917_156_S103", "P11917_159_S106", 
                                                   "P11917_162_S109", "P11917_165_S112"))

Instar5M_coldata <- as.data.frame(colData[c(2, 5, 8, 11, 14, 17, 20, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 61, 64), ])

Instar5F_countdata <- subset(countData, select = c("P11917_103_S51", "P11917_106_S54", "P11917_109_S57", "P11917_112_S60", "P11917_115_S63", 
                                                   "P11917_118_S66", "P11917_121_S69", "P11917_124_S71", "P11917_127_S74", "P11917_130_S77", 
                                                   "P11917_133_S80", "P11917_136_S83", "P11917_139_S86", "P11917_142_S89", "P11917_145_S92", 
                                                   "P11917_148_S95", "P11917_151_S98", "P11917_154_S101", "P11917_157_S104", "P11917_160_S107", 
                                                   "P11917_163_S110", "P11917_166_S113"))
Instar5F_coldata <- as.data.frame(colData[c(3, 6, 9, 12, 15, 18, 21, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65), ])

Instar3_coldata <- Instar3_coldata[ , 1:2]
Instar5M_coldata <- Instar5M_coldata[ , 1:2]
Instar5F_coldata <- Instar5F_coldata[ , 1:2]

#Pre-filtering I: Only keep genes with >1 samples with non-zero read count ####    >>> this is now done above, for all samples (as opposed to doing it for each the individual stage)
#Instar_countdata <- Instar_countdata[rowSums(Instar_countdata > 0) >= 2 , ]
#Pupa_countdata <- Pupa_countdata[rowSums(Pupa_countdata > 0) >= 2 , ]
#Adult_countdata <- Adult_countdata[rowSums(Adult_countdata > 0) >= 2 , ]

#Pre-filtering II: for each developmental stage, remove genes with zero counts in all samples ####  
Instar3_countdata <- Instar3_countdata[rowSums(Instar3_countdata > 0) >= 1 , ]
nrow(Instar3_countdata)
Instar5M_countdata <- Instar5M_countdata[rowSums(Instar5M_countdata > 0) >= 1 , ]
nrow(Instar5M_countdata)
Instar5F_countdata <- Instar5F_countdata[rowSums(Instar5F_countdata > 0) >= 1 , ]
nrow(Instar5F_countdata)


## Add 1 count to every gene/sample (multiple zero-count samples take DESEQ off track)
Instar3_countdata <- Instar3_countdata + 1
head(Instar3_countdata)
Instar5M_countdata <- Instar5M_countdata + 1
head(Instar5M_countdata)
Instar5F_countdata <- Instar5F_countdata + 1
head(Instar5F_countdata)


#Check all sample IDs in colData are also in CountData and match their orders ####
all(rownames(colData) %in% colnames(countData))
countData <- countData[, rownames(colData)]
all(rownames(colData) == colnames(countData))

all(rownames(Instar3_coldata) %in% colnames(Instar3_countdata))
Instar3_countdata <- Instar3_countdata[, rownames(Instar3_coldata)]
all(rownames(Instar3_coldata) == colnames(Instar3_countdata))

all(rownames(Instar5M_coldata) %in% colnames(Instar5M_countdata))
Instar5M_countdata <- Instar5M_countdata[, rownames(Instar5M_coldata)]
all(rownames(Instar5M_coldata) == colnames(Instar5M_countdata))

all(rownames(Instar5F_coldata) %in% colnames(Instar5F_countdata))
Instar5F_countdata <- Instar5F_countdata[, rownames(Instar5F_coldata)]
all(rownames(Instar5F_coldata) == colnames(Instar5F_countdata))




#Create a DESeqDataSet from count matrix and labels (no confounding effects removed at this stage) ####
ddsInstar3 <- DESeqDataSetFromMatrix(countData = Instar3_countdata,
                                    colData = Instar3_coldata, design = ~Population)
ddsInstar3

ddsInstar5M <- DESeqDataSetFromMatrix(countData = Instar5M_countdata,
                                  colData = Instar5M_coldata, design = ~Population)
ddsInstar5M

ddsInstar5F <- DESeqDataSetFromMatrix(countData = Instar5F_countdata,
                                   colData = Instar5F_coldata, design = ~Population)
ddsInstar5F

#Pre-filtering IIIs: Only keep rows with baseMean >2 (usually this woulf be baseMean >1 but we added one count to every gene/sample)
ddsInstar3 <- ddsInstar3[rowMeans(counts(ddsInstar3)) > 2, ]
ddsInstar3
ddsInstar5M <- ddsInstar5M[rowMeans(counts(ddsInstar5M)) > 2, ]
ddsInstar5M
ddsInstar5F <- ddsInstar5F[rowMeans(counts(ddsInstar5F)) > 2, ]
ddsInstar5F











######## MODEL BATCH EFFECTS: Remove hidden batch effects using the sva package

#ddsInstar_sva <- ddsInstar
#ddsInstar_sva <- estimateSizeFactors(ddsInstar_sva)                   # normalization
#dat_Instar  <- counts(ddsInstar_sva, normalized = TRUE)
#idx_Instar  <- rowMeans(dat_Instar) > 2                                 # filter
#dat_Instar  <- dat_Instar[idx_Instar, ]

#ddsPupa_sva <- ddsPupa
#ddsPupa_sva <- estimateSizeFactors(ddsPupa_sva)                   
#dat_Pupa <- counts(ddsPupa_sva, normalized = TRUE)
#idx_Pupa  <- rowMeans(dat_Pupa) > 2                                 
#dat_Pupa  <- dat_Pupa[idx_Pupa, ]

#ddsAdult_sva <- ddsAdult
#ddsAdult_sva <- estimateSizeFactors(ddsAdult_sva)                  
#dat_Adult <- counts(ddsAdult_sva, normalized = TRUE)
#idx_Adult  <- rowMeans(dat_Adult) > 2                                 
#dat_Adult  <- dat_Adult[idx_Adult, ]

# declare model with variable of interest
#mod_Instar  <- model.matrix(~ Sex, colData(ddsInstar_sva))  # when looking for DGE across sexes
#mod_Pupa  <- model.matrix(~ Sex, colData(ddsPupa_sva))
#mod_Adult  <- model.matrix(~ Sex, colData(ddsAdult_sva))

#declare known adjustment variables, or let sva discover batch effects (creates new surrogate variables)
#mod0_Instar <- model.matrix(~ 1, colData(ddsInstar_sva))
#mod0_Pupa <- model.matrix(~ 1, colData(ddsPupa_sva))
#mod0_Adult <- model.matrix(~ 1, colData(ddsAdult_sva))

# estimate number of surrogate variables
#svseq <- svaseq(dat, mod, mod0, n.sv = 3)  ## >>  USE for 03_SDP_LDP (ONLY)
#svseq <- svaseq(dat, mod, mod0, n.sv = 1)  ## >>  USE for 12_LDAF_LDAM (ONLY)
#svseq_Instar <- svaseq(dat_Instar, mod_Instar, mod0_Instar)    ##   If  the sva function  is  called  without  the n.sv argument  specified, the  number  
##   of factors will be estimated for you
#svseq_Instar$sv

#svseq_Pupa <- svaseq(dat_Pupa, mod_Pupa, mod0_Pupa)
#svseq_Pupa$sv

#svseq_Adult <- svaseq(dat_Adult, mod_Adult, mod0_Adult)
#svseq_Adult <- svaseq(dat_Adult, mod_Adult, mod0_Adult, n.sv = 1)  ## >>  USE for 12_LDAF_LDAM (ONLY)
#svseq_Adult$sv

## Use sva to remove any effect of surrogate variables
## (adjust according to number of surrogate variables)
#ddssva_Instar_f <- ddsInstar
#ddssva_Instar_f$SV1 <- svseq_Instar$sv[,1]
#ddssva_Instar_f$SV2 <- svseq_Instar$sv[,2]

#design(ddssva_Instar_f) <- ~ SV1 + SV2 + Sex 

#ddssva_Pupa_f <- ddsPupa
#ddssva_Pupa_f$SV1 <- svseq_Pupa$sv[,1]
#ddssva_Pupa_f$SV2 <- svseq_Pupa$sv[,2]

#design(ddssva_Pupa_f) <- ~ SV1 + SV2 + Sex 

#ddssva_Adult_f <- ddsAdult
#ddssva_Adult_f$SV1 <- svseq_Adult$sv[,1]
#ddssva_Adult_f$SV2 <- svseq_Adult$sv[,2]

#design(ddssva_Adult_f) <- ~ SV1 + SV2 + Sex 
#design(ddssva_Adult_f) <- ~ SV1 + Sex 

             

               











######### Run the default analysis for DESeq2 ####



# Without batch effect correction

ddsInstar3_Pop <- DESeq(ddsInstar3)
res_Instar3_Pop <- results(ddsInstar3_Pop, alpha=.05)
res_Instar3_Pop

ddsInstar5M_Pop <- DESeq(ddsInstar5M)
res_Instar5M_Pop <- results(ddsInstar5M_Pop, alpha=.05)
res_Instar5M_Pop

ddsInstar5F_Pop <- DESeq(ddsInstar5F)
res_Instar5F_Pop <- results(ddsInstar5F_Pop, alpha=.05)
res_Instar5F_Pop



#  With sva correction for batch effects
ddsInstar <- DESeq(ddssva_Instar_f)
res_Instar <- results(ddsInstar, alpha=.05)
res_Instar

ddsPupa <- DESeq(ddssva_Pupa_f)
res_Pupa <- results(ddsPupa, alpha=.05)
res_Pupa

ddsAdult <- DESeq(ddssva_Adult_f)
res_Adult <- results(ddsAdult, alpha=.05)
res_Adult



#Results ####

##Independent hypothesis weighting
#resInstarIHW <- results(ddsInstar, filterFun = ihw, alpha = 0.05)
#resPupaIHW <- results(ddsPupa, filterFun = ihw, alpha = 0.05)
#resAdultIHW <- results(ddsAdult, filterFun = ihw, alpha = 0.05)

##Sort by adjusted p-value and display ####
#(resInstarOrdered <- resInstarIHW[order(resInstarIHW$padj), ])
#(resPupaOrdered <- resPupaIHW[order(resPupaIHW$padj), ])
#(resAdultOrdered <- resAdultIHW[order(resAdultIHW$padj), ])

##summary(resInstar)
#summary(resInstarIHW)
##summary(resPupa)
#summary(resPupaIHW)
##summary(resAdult)
#summary(resAdultIHW)



# Filter p-adj (FDR) < 0.05 ####
#Filtered_P0.05_Instar <- subset(resInstarIHW, padj<0.05)
#Filtered_P0.05_Pupa <- subset(resPupaIHW, padj<0.05)
#Filtered_P0.05_Adult <- subset(resAdultIHW, padj<0.05)


## no batch effects correction

Filtered_P0.05_Instar3_Pop <- subset(res_Instar3_Pop, padj<0.05)
Filtered_P0.05_Instar5M_Pop <- subset(res_Instar5M_Pop, padj<0.05)
Filtered_P0.05_Instar5F_Pop <- subset(res_Instar5F_Pop, padj<0.05)

summary(Filtered_P0.05_Instar3_Pop)
summary(Filtered_P0.05_Instar5M_Pop)
summary(Filtered_P0.05_Instar5F_Pop)


## when correting for batch effects
#Filtered_P0.05_Instar <- subset(res_Instar, padj<0.05)
#Filtered_P0.05_Pupa <- subset(res_Pupa, padj<0.05)
#Filtered_P0.05_Adult <- subset(res_Adult, padj<0.05)

#summary(Filtered_P0.05_Instar)
#summary(Filtered_P0.05_Pupa)
#summary(Filtered_P0.05_Adult)


# Filter logFC > |1.0|) ####
Filtered_P0.05_LFC_1.5_Instar3_Pop <- subset(Filtered_P0.05_Instar3_Pop, log2FoldChange >1 | log2FoldChange < -1)
Filtered_P0.05_LFC_1.5_Instar5M_Pop <- subset(Filtered_P0.05_Instar5M_Pop, log2FoldChange >1 | log2FoldChange < -1)
Filtered_P0.05_LFC_1.5_Instar5F_Pop <- subset(Filtered_P0.05_Instar5F_Pop, log2FoldChange >1 | log2FoldChange < -1)

summary(Filtered_P0.05_LFC_1.5_Instar3_Pop)
summary(Filtered_P0.05_LFC_1.5_Instar5M_Pop)
summary(Filtered_P0.05_LFC_1.5_Instar5F_Pop)


#Filter baseMean > 10 ####
#Filtered_P0.05_LFC_1.5_base_Mean_10_Instar <- subset(Filtered_P0.05_LFC_1.5_Instar, baseMean>10)
#Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa <- subset(Filtered_P0.05_LFC_1.5_Pupa, baseMean>10)
#Filtered_P0.05_LFC_1.5_base_Mean_10_Adult <- subset(Filtered_P0.05_LFC_1.5_Adult, baseMean>10)

#summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar)
#summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa)
#summary(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult)




# Select Male/Female biased genes ####
#Filtered_P0.05_MBG_Instar <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar, log2FoldChange >1)
#Filtered_P0.05_MBG_Pupa <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa , log2FoldChange > 1)
#Filtered_P0.05_MBG_Adult <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, log2FoldChange > 1)
#summary(Filtered_P0.05_MBG_Instar)
#summary(Filtered_P0.05_MBG_Pupa)
#summary(Filtered_P0.05_MBG_Adult)

#Filtered_P0.05_FBG_Instar <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Instar, log2FoldChange < -1)
#Filtered_P0.05_FBG_Pupa <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Pupa, log2FoldChange < -1)
#Filtered_P0.05_FBG_Adult <- subset(Filtered_P0.05_LFC_1.5_base_Mean_10_Adult, log2FoldChange < -1)
#summary(Filtered_P0.05_FBG_Instar)
#summary(Filtered_P0.05_FBG_Pupa)
#summary(Filtered_P0.05_FBG_Adult)

# Filter out genes without expression
#expressed_instar <- subset(resInstarIHW, baseMean >0)
#expressed_pupa <- subset(resPupaIHW, baseMean >0)
#expressed_adult <- subset(resAdultIHW, baseMean >0)

expressed_instar <- subset(Filtered_P0.05_LFC_1.5_Instar3_0, baseMean >0)
expressed_pupa <- subset(Filtered_P0.05_LFC_1.5_Instar5M_0, baseMean >0)
expressed_adult <- subset(Filtered_P0.05_LFC_1.5_Instar5F_0, baseMean >0)

summary(Filtered_P0.05_LFC_1.5_Instar3_0)
summary(Filtered_P0.05_LFC_1.5_Instar5M_0)
summary(Filtered_P0.05_LFC_1.5_Instar5F_0)

# Non-biased genes
Nonbiased_genes_Instar <- subset(res_Instar, padj>0.05 | baseMean<10 | (log2FoldChange <1 & log2FoldChange > -1))
Nonbiased_genes_Pupa <- subset(res_Pupa, padj>0.05 | baseMean<10 | (log2FoldChange <1 & log2FoldChange > -1))
Nonbiased_genes_Adult <- subset(res_Adult, padj>0.05 | baseMean<10 | (log2FoldChange <1 & log2FoldChange > -1))
summary(Nonbiased_genes_Instar)
summary(Nonbiased_genes_Pupa)
summary(Nonbiased_genes_Adult)

Nonbiased_genes_Instar_f <- subset(res_Instar, padj>0.05 | baseMean<10 | log2FoldChange > -1)
Nonbiased_genes_Pupa_f <- subset(res_Pupa, padj>0.05 | baseMean<10 | log2FoldChange > -1)
Nonbiased_genes_Adult_f <- subset(res_Adult, padj>0.05 | baseMean<10 | log2FoldChange > -1)
summary(Nonbiased_genes_Instar_f)
summary(Nonbiased_genes_Pupa_f)
summary(Nonbiased_genes_Adult_f)

Nonbiased_genes_Instar_m <- subset(res_Instar, padj>0.05 | baseMean<10 | log2FoldChange < 1)
Nonbiased_genes_Pupa_m <- subset(res_Pupa, padj>0.05 | baseMean<10 | log2FoldChange < 1)
Nonbiased_genes_Adult_m <- subset(res_Adult, padj>0.05 | baseMean<10 | log2FoldChange < 1)
summary(Nonbiased_genes_Instar_m)
summary(Nonbiased_genes_Pupa_m)
summary(Nonbiased_genes_Adult_m)







###### Write results to file ####


######### standard filtering (genes with p-adj < 0.05, |logFC| > 1, base_Mean > 10)

setwd(OUT_FOLD_STDF)

#unfiltered
write.table(res_Instar3_Pop,
          file = "res_Instar3_Pop.txt", sep = "\t", quote = FALSE)
write.table(res_Instar5M_Pop, 
          file = "res_Instar5M_Pop.txt", sep = "\t", quote = FALSE)
write.table(res_Instar5F_Pop, 
          file = "res_Instar5F_Pop.txt", sep = "\t", quote = FALSE)

#filtered p<0.05
write.table(Filtered_P0.05_Instar3_Pop, 
          file = "Filtered_P0.05_Instar3_Pop.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_Instar5M_Pop, 
          file = "Filtered_P0.05_Instar5M_Pop.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_Instar5F_Pop, 
          file = "Filtered_P0.05_Instar5F_Pop.txt", sep = "\t", quote = FALSE)

#filtered LFC |1|
write.table(Filtered_P0.05_LFC_1.5_Instar3_Pop, 
          file = "Filtered_P0.05_LFC_1.5_Instar3_Pop.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_LFC_1.5_Instar5M_Pop, 
          file = "Filtered_P0.05_LFC_1.5_Instar5M_Pop.txt", sep = "\t", quote = FALSE)
write.table(Filtered_P0.05_LFC_1.5_Instar5F_Pop, 
          file = "Filtered_P0.05_LFC_1.5_Instar5F_Pop.txt", sep = "\t", quote = FALSE)


######################

#All LFC
write.table(expressed_instar, file = "LFC_Instar.txt", sep = "\t", quote = FALSE)
write.table(expressed_pupa, file = "LFC_Pupa.txt", sep = "\t", quote = FALSE)
write.table(expressed_adult, file = "LFC_Adult.txt", sep = "\t", quote = FALSE)


# Non-biased genes
write.table(Nonbiased_genes_Instar, file = "Nonbiased_genes_Instar.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Pupa, file = "Nonbiased_genes_Pupa.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Adult, file = "Nonbiased_genes_Adult.txt", sep = "\t", quote = FALSE)

write.table(Nonbiased_genes_Instar_f, file = "Nonbiased_genes_Instar_f.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Pupa_f, file = "Nonbiased_genes_Pupa_f.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Adult_f, file = "Nonbiased_genes_Adult_f.txt", sep = "\t", quote = FALSE)

write.table(Nonbiased_genes_Instar_m, file = "Nonbiased_genes_Instar_m.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Pupa_m, file = "Nonbiased_genes_Pupa_m.txt", sep = "\t", quote = FALSE)
write.table(Nonbiased_genes_Adult_m, file = "Nonbiased_genes_Adult_m.txt", sep = "\t", quote = FALSE)




######### minimal filtering (all genes with p-adj < 0.05)

#setwd(OUT_FOLD_MINF)



# Select Male/Female biased genes ####
#FilteredMIN_P0.05_MBG_Instar <- subset(Filtered_P0.05_Instar, log2FoldChange >0)
#FilteredMIN_P0.05_MBG_Pupa <- subset(Filtered_P0.05_Pupa, log2FoldChange > 0)
#FilteredMIN_P0.05_MBG_Adult <- subset(Filtered_P0.05_Adult, log2FoldChange > 0)
#summary(FilteredMIN_P0.05_MBG_Instar)
#summary(FilteredMIN_P0.05_MBG_Pupa)
#summary(FilteredMIN_P0.05_MBG_Adult)

#FilteredMIN_P0.05_FBG_Instar <- subset(Filtered_P0.05_Instar, log2FoldChange < 0)
#FilteredMIN_P0.05_FBG_Pupa <- subset(Filtered_P0.05_Pupa, log2FoldChange < 0)
#FilteredMIN_P0.05_FBG_Adult <- subset(Filtered_P0.05_Adult, log2FoldChange < 0)
#summary(FilteredMIN_P0.05_FBG_Instar)
#summary(FilteredMIN_P0.05_FBG_Pupa)
#summary(FilteredMIN_P0.05_FBG_Adult)

# Non-biased genes >>>> NOT WORKING AS CODE BELOW DOES NOT GREP GENES WITH padj=NA
#Nonbiased_genes_Instar_MIN <- subset(res_Instar, padj>0.05)
#Nonbiased_genes_Pupa_MIN <- subset(res_Pupa, padj>0.05)
#Nonbiased_genes_Adult_MIN <- subset(res_Adult, padj>0.05)
#summary(Nonbiased_genes_Instar_MIN)
#summary(Nonbiased_genes_Pupa_MIN)
#summary(Nonbiased_genes_Adult_MIN)

#Nonbiased_genes_Instar_MIN_f <- subset(res_Instar, padj>0.05 | log2FoldChange > 0)
#Nonbiased_genes_Pupa_MIN_f <- subset(res_Pupa, padj>0.05 | log2FoldChange > 0)
#Nonbiased_genes_Adult_MIN_f <- subset(res_Adult, padj>0.05 | log2FoldChange > 0)
#summary(Nonbiased_genes_Instar_MIN_f)
#summary(Nonbiased_genes_Pupa_MIN_f)
#summary(Nonbiased_genes_Adult_MIN_f)

#Nonbiased_genes_Instar_MIN_m <- subset(res_Instar, padj>0.05 | log2FoldChange < 0)
#Nonbiased_genes_Pupa_MIN_m <- subset(res_Pupa, padj>0.05 | log2FoldChange < 0)
#Nonbiased_genes_Adult_MIN_m <- subset(res_Adult, padj>0.05 | log2FoldChange < 0)
#summary(Nonbiased_genes_Instar_MIN_m)
#summary(Nonbiased_genes_Pupa_MIN_m)
#summary(Nonbiased_genes_Adult_MIN_m)



#write.table(Filtered_P0.05_Instar,
#            file = "Filtered_P0.05_Instar.txt", sep = "\t", quote = FALSE)
#write.table(Filtered_P0.05_Pupa, 
#            file = "Filtered_P0.05_Pupa.txt", sep = "\t", quote = FALSE)
#write.table(Filtered_P0.05_Adult, 
#            file = "Filtered_P0.05_Adult.txt", sep = "\t", quote = FALSE)

#write.table(FilteredMIN_P0.05_MBG_Instar, 
#            file = "MBG_Instar.txt", sep = "\t", quote = FALSE)
#write.table(FilteredMIN_P0.05_MBG_Pupa, 
#            file = "MBG_Pupa.txt", sep = "\t", quote = FALSE)
#write.table(FilteredMIN_P0.05_MBG_Adult, 
#            file = "MBG_Adult.txt", sep = "\t", quote = FALSE)
#
#write.table(FilteredMIN_P0.05_FBG_Instar, 
#            file = "FBG_Instar.txt", sep = "\t", quote = FALSE)
#write.table(FilteredMIN_P0.05_FBG_Pupa, 
#            file = "FBG_Pupa.txt", sep = "\t", quote = FALSE)
#write.table(FilteredMIN_P0.05_FBG_Adult, 
#            file = "FBG_Adult.txt", sep = "\t", quote = FALSE)

#All LFC
#write.table(expressed_instar, file = "LFC_Instar.txt", sep = "\t", quote = FALSE)
#write.table(expressed_pupa, file = "LFC_Pupa.txt", sep = "\t", quote = FALSE)
#write.table(expressed_adult, file = "LFC_Adult.txt", sep = "\t", quote = FALSE)

# Non-biased genes
#write.table(Nonbiased_genes_Instar_MIN, file = "Nonbiased_genes_Instar.txt", sep = "\t", quote = FALSE)
#write.table(Nonbiased_genes_Pupa_MIN, file = "Nonbiased_genes_Pupa.txt", sep = "\t", quote = FALSE)
#write.table(Nonbiased_genes_Adult_MIN, file = "Nonbiased_genes_Adult.txt", sep = "\t", quote = FALSE)

#write.table(Nonbiased_genes_Instar_MIN_f, file = "Nonbiased_genes_Instar_f.txt", sep = "\t", quote = FALSE)
#write.table(Nonbiased_genes_Pupa_MIN_f, file = "Nonbiased_genes_Pupa_f.txt", sep = "\t", quote = FALSE)
#write.table(Nonbiased_genes_Adult_MIN_f, file = "Nonbiased_genes_Adult_f.txt", sep = "\t", quote = FALSE)

#write.table(Nonbiased_genes_Instar_MIN_m, file = "Nonbiased_genes_Instar_m.txt", sep = "\t", quote = FALSE)
#write.table(Nonbiased_genes_Pupa_MIN_m, file = "Nonbiased_genes_Pupa_m.txt", sep = "\t", quote = FALSE)
#write.table(Nonbiased_genes_Adult_MIN_m, file = "Nonbiased_genes_Adult_m.txt", sep = "\t", quote = FALSE)

######### 


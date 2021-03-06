
###########
# Created by Karin Näsvall, Uppsala universitet
# Project RNAseq_microbiome_hostplant_project
# 2019-12-13
###########

#This script makes a topGO analysis of differetially expressed genes, used in the Leptidea_Hostplant-project 2019
#obs!! 
#rm(list=ls())  #run it only in console before running the script, should not be used in Rmarkdown


library(topGO)
library(AnnotationDbi)
library(GO.db)
library(org.Dm.eg.db)

Sys.Date()
sessionInfo()


##################################
#get the data
#create gene list from list created with Stringtie where the gene alias in orthologs to D. melanogaster is included. 
#create a background set of genes from the specific dataset filtered to include the genes with basemean > 10
#select genes that are differentially expresssed from the different DE-datasets 

########
# OBS!!! Change "instar5Mspan" in the whole document to the name for the actual dataset
########

#background genes for the dataset, genes are filtered to over 10 reads basemean
background_input_file <- "../../Results/4_DEseq2/06_Spa_Hostplant_effect/Filtered_base_Mean_10_Instar5M_Spa_Host.txt" #CHANGE THIS!
#candidate genes with significant different geneexpression, the geneset of interest
candidate_input_file <- "../../Results/4_DEseq2/06_Spa_Hostplant_effect/Filtered_base_Mean_10_LFC_1.0_P0.05_Instar5M_Spa_Host.txt" #CHANGE THIS!


#create gene list with geneID this can be reused for all datasets #Shouls be CHANGED if workning with other organisms or databases
geneList_dmel=read.table("gene_list_sin_dmel_stringtie.txt", header = T)
#get rid of the "-" (unannotated genes) 
geneList_dmel_filtered <- subset(geneList_dmel, !geneList_dmel$GeneName== "-")
length(geneList_dmel_filtered$GeneID)


#read input background genes
bg_instar5Mspan <- read.table(background_input_file, header = T)
bg_instar5Mspan$GeneID <- rownames(bg_instar5Mspan)
#get nr og genes
length(bg_instar5Mspan$GeneID)

#read input DE genes 
candidateGenes_instar5Mspan <- read.table(candidate_input_file, header = T)
candidateGenes_instar5Mspan$GeneID <- rownames(candidateGenes_instar5Mspan)
#get nr og genes
length(candidateGenes_instar5Mspan$GeneID)


#get the D mel gene names
bg_instar5Mspan_Dmel <- merge(bg_instar5Mspan, geneList_dmel_filtered, by.x="GeneID", by.y="GeneID")
candidateGenes_instar5Mspan_Dmel <- merge(candidateGenes_instar5Mspan, 
                                          geneList_dmel_filtered, 
                                          by.x="GeneID", 
                                          by.y="GeneID")
#manually curating GeneName list
write(as.character(candidateGenes_instar5Mspan_Dmel$GeneName), file="candidateGenes_instar5Mspan_Dmel_prel.txt")
#get new list for manually curating the genenames (ex amyamy35amy2 = amy)
candidateGenes_instar5Mspan_Dmel_new <- read.table("candidateGenes_instar5Mspan_Dmel.txt", col.names = "GeneName")
candidateGenes_instar5Mspan_Dmel_new_list <- as.character(candidateGenes_instar5Mspan_Dmel_new$GeneName)
#length(candidateGenes_instar5Mspan_Dmel$GeneID)
#the gene names have to be in a vector as characters
candidateGenes_instar5Mspan_old_list <- as.character(candidateGenes_instar5Mspan_Dmel$GeneName)
length(candidateGenes_instar5Mspan_old_list)

#get the names that are different between the two lists, and number of genenames changed
candidateGenes_instar5Mspan_Dmel_new_diff <- candidateGenes_instar5Mspan_Dmel_new_list[!(candidateGenes_instar5Mspan_Dmel_new_list %in% candidateGenes_instar5Mspan_old_list)]
length(candidateGenes_instar5Mspan_Dmel_new_diff)

######################

bg_instar5Mspan_Dmel_list <- as.character(bg_instar5Mspan_Dmel$GeneName)
length(bg_instar5Mspan_Dmel_list)

bg_instar5Mspan_Dmel_new_list <- c(bg_instar5Mspan_Dmel_list,candidateGenes_instar5Mspan_Dmel_new_diff)
length(bg_instar5Mspan_Dmel_new_list)

###########################################
#GO analysis
# GeneName to GO mapping, get nr of GO terms per gene in the candadate genelist
GOtoalias_instar5Mspan_BP <- annFUN.org("BP", 
                                        feasibleGenes = candidateGenes_instar5Mspan_Dmel_new_list, 
                                        mapping = "org.Dm.eg.db", 
                                        ID = "Alias")
GeneIDtoGO_instar5Mspan_BP <- inverseList(GOtoalias_instar5Mspan_BP) # only 2 genes with GO annotations
GeneIDtoGO_instar5Mspan_BP

GOtoalias_instar5Mspan_MF <- annFUN.org("MF", 
                                        feasibleGenes = candidateGenes_instar5Mspan_Dmel_new_list, 
                                        mapping = "org.Dm.eg.db", 
                                        ID = "Alias")
GeneIDtoGO_instar5Mspan_MF <- inverseList(GOtoalias_instar5Mspan_MF) # only 2 genes with GO annotations
GeneIDtoGO_instar5Mspan_MF

GOtoalias_instar5Mspan_CC <- annFUN.org("CC", 
                                        feasibleGenes = candidateGenes_instar5Mspan_Dmel_new_list, 
                                        mapping = "org.Dm.eg.db", 
                                        ID = "Alias")
GeneIDtoGO_instar5Mspan_CC <- inverseList(GOtoalias_instar5Mspan_CC) # only 2 genes with GO annotations
GeneIDtoGO_instar5Mspan_CC

#convert these list into a dataframe to be able to merge it later in the result folder, add a column so that gene-aliases can be a separate column
GeneIDtoGO_instar5Mspan_BP.df <- do.call(rbind, lapply(GeneIDtoGO_instar5Mspan_BP, as.data.frame))
colnames(GeneIDtoGO_instar5Mspan_BP.df) <- c("BP")
GeneIDtoGO_instar5Mspan_BP.df$GeneName <- rownames(GeneIDtoGO_instar5Mspan_BP.df)

GeneIDtoGO_instar5Mspan_MF.df <- do.call(rbind, lapply(GeneIDtoGO_instar5Mspan_MF, as.data.frame))
colnames(GeneIDtoGO_instar5Mspan_MF.df) <- c("MF")
GeneIDtoGO_instar5Mspan_MF.df$GeneName <- rownames(GeneIDtoGO_instar5Mspan_MF.df)

GeneIDtoGO_instar5Mspan_CC.df <- do.call(rbind, lapply(GeneIDtoGO_instar5Mspan_CC, as.data.frame))
colnames(GeneIDtoGO_instar5Mspan_CC.df) <- c("CC")
GeneIDtoGO_instar5Mspan_CC.df$GeneName <- rownames(GeneIDtoGO_instar5Mspan_CC.df)

#############
#creating an topGOdata-object
#first give the candidate genes a "score" Listed as 1 and the other are listed as 0
geneList_instar5Mspan <- factor(as.integer(bg_instar5Mspan_Dmel_new_list %in% candidateGenes_instar5Mspan_Dmel_new_list))
names(geneList_instar5Mspan) <- bg_instar5Mspan_Dmel_new_list
str(geneList_instar5Mspan)

sampleGOdata_instar5Mspan_BP <- new(Class = "topGOdata", 
                                    description="instar5Mspan_hostplant", 
                                    ontology="BP", 
                                    allGenes=geneList_instar5Mspan, 
                                    annot=annFUN.org, 
                                    mapping = "org.Dm.eg.db", 
                                    ID = "Alias")
#get a short summary
sampleGOdata_instar5Mspan_BP

sampleGOdata_instar5Mspan_MF <- new(Class = "topGOdata", 
                                    description="instar5Mspan_hostplant", 
                                    ontology="MF", 
                                    allGenes=geneList_instar5Mspan, 
                                    annot=annFUN.org, 
                                    mapping = "org.Dm.eg.db", 
                                    ID = "Alias")
#get a short summary
sampleGOdata_instar5Mspan_MF

sampleGOdata_instar5Mspan_CC <- new(Class = "topGOdata", 
                                    description="instar5Mspan_hostplant", 
                                    ontology="CC", 
                                    allGenes=geneList_instar5Mspan, 
                                    annot=annFUN.org, 
                                    mapping = "org.Dm.eg.db", 
                                    ID = "Alias")
#get a short summary
sampleGOdata_instar5Mspan_CC




###########################################################
#run a significance test 

######## BP
#The elim and weight algorithms were introduced in Alexa et al. (2006). The default algorithm used by the topGO package is a mixture between #the elim and the weight algorithms and it will be referred as weight01.
weight_fisher_result_instar5Mspan_BP=runTest(sampleGOdata_instar5Mspan_BP, 
                                             algorithm='weight01', 
                                             statistic='fisher')
#get the result summary
weight_fisher_result_instar5Mspan_BP
#get the p-values with score
#weight_fisher_result_pvalues_instar5Mspan <- score(weight_fisher_result_instar5Mspan)
#more results
#geneData(weight_fisher_result_instar5Mspan_BP)
#get results in table
results_table_instar5Mspan_BP <- GenTable(object = sampleGOdata_instar5Mspan_BP, 
                                          weight=weight_fisher_result_instar5Mspan_BP, 
                                          orderBy="weight", 
                                          topNodes=20, 
                                          numChar=100)
#results_table_instar5Mspan_BP

#get only GOterms sign enriched , sorted by GO term
#sort(subset(results_table_instar5Mspan, results_table_instar5Mspan$weight <= 0.05)[ ,1])
#get only significantly enriched GOterm with result table
results_table_instar5Mspan_sign_BP <- subset(results_table_instar5Mspan_BP, results_table_instar5Mspan_BP$weight <= 0.05)
results_table_instar5Mspan_sign_BP

#checking the distribution of annotations
#results_table_instar5Mspan_long <- GenTable(object = sampleGOdata_instar5Mspan_BP, weight=weight_fisher_result_instar5Mspan, orderBy="weight", topNodes=1000)
#hist(results_table_instar5Mspan_long$Annotated, breaks = 500, xlim = c(0,100))

####### MF
weight_fisher_result_instar5Mspan_MF=runTest(sampleGOdata_instar5Mspan_MF, algorithm='weight01', statistic='fisher')
weight_fisher_result_instar5Mspan_MF
results_table_instar5Mspan_MF <- GenTable(object = sampleGOdata_instar5Mspan_MF, 
                                          weight=weight_fisher_result_instar5Mspan_MF, 
                                          orderBy="weight", 
                                          topNodes=20, 
                                          numChar=100)

#get only significantly enriched GOterm with result table
results_table_instar5Mspan_sign_MF <- subset(results_table_instar5Mspan_MF, results_table_instar5Mspan_MF$weight <= 0.05)
results_table_instar5Mspan_sign_MF

####### CC
weight_fisher_result_instar5Mspan_CC=runTest(sampleGOdata_instar5Mspan_CC, algorithm='weight01', statistic='fisher')
weight_fisher_result_instar5Mspan_CC
results_table_instar5Mspan_CC <- GenTable(object = sampleGOdata_instar5Mspan_CC, 
                                          weight=weight_fisher_result_instar5Mspan_CC, 
                                          orderBy="weight", 
                                          topNodes=20, 
                                          numChar=100)

#get only significantly enriched GOterm with result table
results_table_instar5Mspan_sign_CC <- subset(results_table_instar5Mspan_CC, results_table_instar5Mspan_CC$weight <= 0.05)
results_table_instar5Mspan_sign_CC

##########################################
#collecting the results in one table
#GeneID and GOterms assoc with those genes, GeneIDtoGO_instar5Mspan_BP etc, 
#GOterms and p-values (table), results_table_instar5Mspan_sign_BP

#get genenames from the org....db, can be resused in the different datasets
GeneListNames_db <- select(org.Dm.eg.db, keys=keys(org.Dm.eg.db), columns=c("GENENAME", "ALIAS"))

###### BP  #####################
result_table_GO_instar5Mspan_BP <- merge(results_table_instar5Mspan_sign_BP, 
                                         GeneIDtoGO_instar5Mspan_BP.df, 
                                         by.x="GO.ID", 
                                         by.y="BP")
#remove the extra suffix ork1.1 etc
result_table_GO_instar5Mspan_BP$GeneName=gsub(pattern = "\\.[0-99]*$",  
                                              replacement = "", 
                                              x = result_table_GO_instar5Mspan_BP$GeneName, 
                                              fixed = FALSE)
#result_table_GO_instar5Mspan_BP

#merge df to one list with GO.ID, GOterms, gene-alias and gene names, and p-value for enrichment-test
result_table_GO_instar5Mspan_BP_temp <- merge(result_table_GO_instar5Mspan_BP, 
                                              GeneListNames_db, 
                                              by.x="GeneName", 
                                              by.y="ALIAS")
#to add a column with GeneID for L.sinapis
result_table_GO_instar5Mspan_BP_genenames <- merge(result_table_GO_instar5Mspan_BP_temp, 
                                                   geneList_dmel_filtered, 
                                                   by.x="GeneName", 
                                                   by.y="GeneName")
result_table_GO_instar5Mspan_BP_genenames


###### MF  ####################
result_table_GO_instar5Mspan_MF <- merge(results_table_instar5Mspan_sign_MF, 
                                         GeneIDtoGO_instar5Mspan_MF.df, 
                                         by.x="GO.ID", 
                                         by.y="MF")
#remove the extra suffix ork1.1 etc
result_table_GO_instar5Mspan_MF$GeneName=gsub(pattern = "\\.[0-99]*$",  
                                              replacement = "", 
                                              x = result_table_GO_instar5Mspan_MF$GeneName, 
                                              fixed = FALSE)
#result_table_GO_instar5Mspan_MF

#merge df to one list with GO.ID, GOterms, gene-alias and gene names, and p-value for enrichment-test
result_table_GO_instar5Mspan_MF_temp <- merge(result_table_GO_instar5Mspan_MF, 
                                              GeneListNames_db, 
                                              by.x="GeneName", 
                                              by.y="ALIAS")
#to add a column with GeneID for L.sinapis
result_table_GO_instar5Mspan_MF_genenames <- merge(result_table_GO_instar5Mspan_MF_temp, 
                                                   geneList_dmel_filtered, 
                                                   by.x="GeneName", 
                                                   by.y="GeneName")
result_table_GO_instar5Mspan_MF_genenames


###### CC  ###################
result_table_GO_instar5Mspan_CC <- merge(results_table_instar5Mspan_sign_CC, 
                                         GeneIDtoGO_instar5Mspan_CC.df, 
                                         by.x="GO.ID", 
                                         by.y="CC")
#remove the extra suffix ork1.1 etc
result_table_GO_instar5Mspan_CC$GeneName=gsub(pattern = "\\.[0-99]*$",  
                                              replacement = "", 
                                              x = result_table_GO_instar5Mspan_CC$GeneName, 
                                              fixed = FALSE)
#result_table_GO_instar5Mspan_CC

#merge df to one list with GO.ID, GOterms, gene-alias and gene names, and p-value for enrichment-test
result_table_GO_instar5Mspan_CC_temp <- merge(result_table_GO_instar5Mspan_CC, 
                                              GeneListNames_db, 
                                              by.x="GeneName", 
                                              by.y="ALIAS")
#to add a column with GeneID for L.sinapis
result_table_GO_instar5Mspan_CC_genenames <- merge(result_table_GO_instar5Mspan_CC_temp, 
                                                   geneList_dmel_filtered, 
                                                   by.x="GeneName", 
                                                   by.y="GeneName")
result_table_GO_instar5Mspan_CC_genenames


#write result to table
write.csv(result_table_GO_instar5Mspan_BP_genenames, file="../../Results/6_topGO/result_GO_instar5Mspan_hostplant_BP.table")
write.csv(result_table_GO_instar5Mspan_MF_genenames, file="../../Results/6_topGO/result_GO_instar5Mspan_hostplant_MF.table")
write.csv(result_table_GO_instar5Mspan_CC_genenames, file="../../Results/6_topGO/result_GO_instar5Mspan_hostplant_CC.table")

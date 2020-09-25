sessionInfo()
library(topGO)
library(org.Dm.eg.db)
library(AnnotationDbi)

#create gene universe from gene list created with Stringtie where the gene names in orthologos to D. melanogaster is included. 
#select genes that are differentially expresssed from the different DE-datasets 

#first create gene universe with geneID 
geneUniv=read.table("gene_list_sin_dmel_stringtie.txt", header = T)

#background genes with only D.mel names
bg_genes<-as.character(geneUniv[,2])
#get rid of the "-" (unannotated genes) 
bg_genes_filtered <- as.character(subset(bg_genes, !bg_genes== "-"))

#DE genes in Cat instar3
candidateGenes <- read.table("../../Results/4_DEseq2/06_Spa_Hostplant_effect/Filtered_base_Mean_10_LFC_1.0_P0.05_Instar3_Spa_Host.txt", header = T)
candidateGenes$GeneID <- rownames(candidateGenes)
length(candidateGenes)

#get the D mel gene names for the candidate genes
canditateGenes_Dmel <- merge(candidateGenes, geneUniv, by.x="GeneID", by.y="GeneID")

candidateGenes_filtered <- as.character(subset(canditateGenes_Dmel$GeneName, !canditateGenes_Dmel$GeneName=="-"))
length(candidateGenes_filtered) #33 genes
length(bg_genes_filtered) #10926
# GeneName to GO mapping needed. 


GOtoalias_dmel <- annFUN.org("BF", feasibleGenes = candidateGenes_filtered, mapping = "org.Dm.eg.db", ID = "Alias")
GeneIDtoGOex <- inverseList(GOtoalias_dmel) # only 6 genes with annotations
GOtoalias_dmel_MF <- annFUN.org("MF", feasibleGenes = candidateGenes_filtered, mapping = "org.Dm.eg.db", ID = "Alias")
GeneIDtoGOex_MF <- inverseList(GOtoalias_dmel_MF) # only 6 genes with annotations

GOtoalias_dmel_CC <- annFUN.org("CC", feasibleGenes = candidateGenes_filtered, mapping = "org.Dm.eg.db", ID = "Alias")
GeneIDtoGOex_CC <- inverseList(GOtoalias_dmel_CC) # only 6 genes with annotations


geneList_ex <- factor(as.integer(bg_genes_filtered %in% candidateGenes_filtered))
names(geneList_ex) <- bg_genes_filtered
str(geneList_ex)

sampleGOdata_ex_BP <- new(Class = "topGOdata", description="test_session", ontology="BP", allGenes=geneList_ex, annot=annFUN.org, mapping = "org.Dm.eg.db", ID = "Alias")
#get a short summary
sampleGOdata_ex_BP
#run a significance test
weight_fisher_result=runTest(sampleGOdata_ex_BP, algorithm='weight01', statistic='fisher')
#get the result
weight_fisher_result
#get the p-values with score
weight_fisher_result_pvalues <- score(weight_fisher_result)
#more results
geneData(weight_fisher_result)
#Write results in table
results_table_ex <- GenTable(object = sampleGOdata_ex_BP, weight=weight_fisher_result, orderBy="weight", topNodes=20)
results_table_ex
sort(results_table_ex[ ,1])
#get only sign enriched GOterms
sort(subset(results_table_ex, results_table_ex$weight <= 0.05)[ ,1])

#checking the distribution of annotations
results_table_ex_long <- GenTable(object = sampleGOdata_ex_BP, weight=weight_fisher_result, orderBy="weight", topNodes=1000)
hist(results_table_ex_long$Annotated, breaks = 500, xlim = c(0,100))
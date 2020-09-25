#r script
#Visualising the data, obs can be used by copy and paste the commands into the r script where the data is loaded and DEseq2 analysis is run.

#plotting the results 
#MA plot (but Volcanoplot is better)
par(mfrow=c(2,2))
plotMA(res_All_Body_contr_plant, main= "All, Lot vs Dor")
plotMA(res_Instar3_Pop_Host, main="Instar3, Lot vs Dor")
plotMA(res_Instar5M_Pop_Host, main="Instar5M, Lot vs Dor")
plotMA(res_Instar5F_Pop_Host, main="Instar5F, Lot vs Dor")
plotMA(res_All_Body_contr, main="All, span vs swe")
plotMA(res_Instar3_contr1_Pop, main="Instar3, span vs swe")
plotMA(res_Instar5M_contr1_Pop, main="Instar5M, span vs swe")
plotMA(res_Instar5F_contr1_Pop, main="Instar5F, span vs swe")

plotCounts(ddsAll_Body, gene=which.min(res_All_Body$padj), intgroup="Hostplant", xlab = "All")
plotCounts(ddsInstar3_Pop_Host, gene=which.min(res_Instar3_Pop_Host$padj), intgroup="Hostplant", xlab = "Instar 3")
plotCounts(ddsInstar5M_Pop_Host, gene=which.min(res_Instar5M_Pop_Host$padj), intgroup="Hostplant", xlab = "Instar 5 Male")
plotCounts(ddsInstar5F_Pop_Host, gene=which.min(res_Instar5F_Pop_Host$padj), intgroup="Hostplant", xlab = "Instar 5 Female")


#Heatmap
install.packages("pheatmap")
library("pheatmap")

#All
All_Body.ntd <- normTransform(ddsAll_Body)
select_All_Body <- order(rowMeans(counts(ddsAll_Body,normalized=TRUE)),
                         decreasing=TRUE)

All_Body.df <- as.data.frame(colData(ddsAll_Body)[,c("Population", "Hostplant")])
pheatmap(assay(All_Body.ntd)[select_All_Body,order(All_Body.ntd$Hostplant)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=All_Body.df)

#cluster columns
All_Body.df <- as.data.frame(colData(ddsAll_Body)[,c("Population", "Hostplant", "Stage")])
pheatmap(assay(All_Body.ntd)[select_All_Body,order(All_Body.ntd$Hostplant)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=All_Body.df)


#Instar5
Instar5_Body.ntd <- normTransform(ddsInstar5_Body)
select_Instar5_Body <- order(rowMeans(counts(ddsInstar5_Body,normalized=TRUE)),
                             decreasing=TRUE)
Instar5_Body.df <- as.data.frame(colData(ddsInstar5_Body)[,c("Population", "Hostplant")])
pheatmap(assay(Instar5_Body.ntd)[select_Instar5_Body,order(Instar5_Body.ntd$Hostplant)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=Instar5_Body.df)

#Instar3
Instar3.ntd <- normTransform(ddsInstar3_Pop_Host)
select_Instar3 <- order(rowMeans(counts(ddsInstar3_Pop_Host,normalized=TRUE)),
                        decreasing=TRUE)
Instar3.df <- as.data.frame(colData(ddsInstar3_Pop_Host)[,c("Population", "Hostplant")])
pheatmap(assay(Instar3.ntd)[select_Instar3,order(Instar3.ntd$Hostplant)], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=Instar3.df)


#PCA
library(ggplot2)
All_Body.vsd <- varianceStabilizingTransformation(ddsAll_Body, blind = T)
plotPCA(All_Body.vsd, intgroup=c("Hostplant","Population")) + theme_minimal() + ggtitle("PCA All")

par(mfrow=c(2,2))
plotPCA(All_Body.vsd, intgroup=c("Hostplant","Population", "Stage")) + theme_minimal()
plotPCA(All_Body.vsd, intgroup=c("Stage")) + theme_minimal() + ggtitle("PCA All Instar 3 vs 5")
plotPCA(All_Body.vsd, intgroup=c("Population")) + theme_minimal() + ggtitle("PCA All Population")
plotPCA(All_Body.vsd, intgroup=c("Hostplant")) + theme_minimal() + ggtitle("PCA All Hostplant")

Instar5_Body.vsd <- varianceStabilizingTransformation(ddsInstar5_Body, blind = T)
pcaInstar5_body <- plotPCA(Instar5_Body.vsd, intgroup=c("Population", "Hostplant")) + ggtitle("PCA Instar 5")
plotPCA(Instar5_Body.vsd, intgroup=c("Hostplant")) + theme_minimal() + ggtitle("PCA Instar 5 Hostplant")
plotPCA(Instar5_Body.vsd, intgroup=c("Population")) + theme_minimal() + ggtitle("PCA Instar 5 Population")
plotPCA(Instar5_Body.vsd, intgroup=c("Sex")) + theme_minimal() + ggtitle("PCA Instar 5 Sex")
plotPCA(Instar5_Body.vsd, intgroup=c("Family")) + theme_minimal() + ggtitle("PCA Instar 5 Family")

Instar3.vsd <- varianceStabilizingTransformation(ddsInstar3_Pop_Host, blind = T)
plotPCA(Instar3.vsd, intgroup=c("Population", "Hostplant")) + theme_minimal() + ggtitle("PCA Instar 3")
plotPCA(Instar3.vsd, intgroup=c("Hostplant")) + theme_minimal() + ggtitle("PCA Instar 3 Hostplant")
plotPCA(Instar3.vsd, intgroup=c("Population")) + theme_minimal() + ggtitle("PCA Instar 3 Population")
plotPCA(Instar3.vsd, intgroup=c("Family")) + theme_minimal() + ggtitle("PCA Instar 3 Family")

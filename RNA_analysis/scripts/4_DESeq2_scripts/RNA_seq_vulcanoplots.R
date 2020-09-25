#Script for Vulcanoplot of DEseq2 data
#20200319 Karin Näsvall
#
install.packages("ggrepel")

library(airway)
library(magrittr)
library('DESeq2')
library(EnhancedVolcano)
library(gridExtra)
getwd()




#get data for swedish and cataloinan separate
res_catIII <- read.table("../Results2/4_DEseq2/06_Spa_Hostplant_effect/res_Instar3_Spa_Host.txt", header = T)
str(res_catIII)

res_catVf <- read.table("../Results2/4_DEseq2/06_Spa_Hostplant_effect/res_Instar5F_Spa_Host.txt", header = T)
str(res_catVf)

res_catVm <- read.table("../Results2/4_DEseq2/06_Spa_Hostplant_effect/res_Instar5M_Spa_Host.txt", header = T)
str(res_catVm)

res_sweIII <- read.table("../Results2/4_DEseq2/05_Swe_Hostplant_effect/res_Instar3_Swe_Host.txt", header = T)
str(res_sweIII)

res_sweVf <- read.table("../Results2/4_DEseq2/05_Swe_Hostplant_effect/res_Instar5F_Swe_Host.txt", header = T)
str(res_sweVf)

res_sweVm <- read.table("../Results2/4_DEseq2/05_Swe_Hostplant_effect/res_Instar5M_Swe_Host.txt", header = T)
str(res_sweVm)

#plots p-value 0.05 = 5e-2)
vplot_catIII<-EnhancedVolcano(res_catIII,
                lab = rownames(res_catIII),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-7.5, 7.5),
                ylim = c(0,7.5),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 5e-2,
                FCcutoff = 1.0,
                gridlines.major = F,
                gridlines.minor = F,
                selectLab = 0,
                title = "CatIII", 
                subtitle = NULL,
                legendVisible = F)
                # legendLabSize = 16, 
                # legendPosition = "bottom", 
                # legend = c("NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))


  #geom_point(size = 3.0)


vplot_catVf<-EnhancedVolcano(res_catVf,
                              lab = rownames(res_catVf),
                              x = 'log2FoldChange',
                              y = 'padj',
                              xlim = c(-7.5, 7.5),
                              ylim = c(0,7.5),
                              ylab = bquote(~-Log[10]~adjusted~italic(P)),
                              pCutoff = 5e-2,
                              FCcutoff = 1.0,
                              gridlines.major = F,
                              gridlines.minor = F,
                              selectLab = 0,
                              title = "CatVf", 
                              subtitle = NULL,
                             legendVisible = F)
                             #  legendLabSize = 16, 
                             #  legendPosition = "bottom", 
                             # legend = c("NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))

vplot_catVm<-EnhancedVolcano(res_catVm,
                              lab = rownames(res_catVm),
                              x = 'log2FoldChange',
                              y = 'padj',
                              xlim = c(-7.5, 7.5),
                              ylim = c(0,7.5),
                              ylab = bquote(~-Log[10]~adjusted~italic(P)),
                              pCutoff = 5e-2,
                              FCcutoff = 1.0,
                              gridlines.major = F,
                              gridlines.minor = F,
                              selectLab = 0,
                              title = "CatVm", 
                              subtitle = NULL, 
                             legendVisible = F)
                             #  legendLabSize = 16, 
                             #  legendPosition = "bottom", 
                             # legend = c("NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))


vplot_sweIII<-EnhancedVolcano(res_sweIII,
                              lab = rownames(res_sweIII),
                              x = 'log2FoldChange',
                              y = 'padj',
                              xlim = c(-7.5, 7.5),
                              ylim = c(0,7.5),
                              ylab = bquote(~-Log[10]~adjusted~italic(P)),
                              pCutoff = 5e-2,
                              FCcutoff = 1.0,
                              gridlines.major = F,
                              gridlines.minor = F,
                              selectLab = 0,
                              title = "SweIII", 
                              subtitle = NULL,
                              legendVisible = F)
#                               legendLabSize = 16, 
#                               legendPosition = "bottom",
#                               legend = c"NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))
# #geom_point(size = 3.0)

vplot_sweVf<-EnhancedVolcano(res_sweVf,
                             lab = rownames(res_sweVf),
                             x = 'log2FoldChange',
                             y = 'padj',
                             xlim = c(-7.5, 7.5),
                             ylim = c(0,7.5),
                             ylab = bquote(~-Log[10]~adjusted~italic(P)),
                             pCutoff = 5e-2,
                             FCcutoff = 1.0,
                             gridlines.major = F,
                             gridlines.minor = F,
                             selectLab = 0,
                             title = "SweVf", 
                             subtitle = NULL, 
                             legendVisible = F)
                             # legendLabSize = 16, 
                             # legendPosition = "bottom", 
                             # legend = c("NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))

vplot_sweVm<-EnhancedVolcano(res_sweVm,
                             lab = rownames(res_sweVm),
                             x = 'log2FoldChange',
                             y = 'padj',
                             xlim = c(-7.5, 7.5),
                             ylim = c(0,7.5),
                             ylab = bquote(~-Log[10]~adjusted~italic(P)),
                             pCutoff = 5e-2,
                             FCcutoff = 1.0,
                             gridlines.major = F,
                             gridlines.minor = F,
                             selectLab = 0,
                             title = "SweVm", 
                             subtitle = NULL, 
                             legendVisible = F)
                             # legendLabSize = 16, 
                             # legendPosition = "bottom", 
                             # legend = c("NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))

#
#get the legend
EnhancedVolcano(res_sweVm,
                lab = rownames(res_sweVm),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-7.5, 7.5),
                ylim = c(0,7.5),
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                pCutoff = 5e-2,
                FCcutoff = 1.0,
                gridlines.major = F,
                gridlines.minor = F,
                selectLab = 0,
                title = "SweVm", 
                subtitle = NULL, 
                legendLabSize = 16,
                legendPosition = "bottom", 
                legend = c("NS", "Log2FC > 1.0", "Padj < 0.05", "Padj < 0.05 & Log2FC > 1.0"))

grid.arrange(vplot_catIII, vplot_sweIII, nrow=1)
grid.arrange(vplot_catVf, vplot_sweVf, nrow=1)
grid.arrange(vplot_catVm, vplot_sweVm, nrow=1)

#, vplot_catVf, vplot_sweVf, vplot_catVm, vplot_sweVm, nrow = 3)

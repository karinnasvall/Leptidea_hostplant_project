###############
#
# @author: Axel Künstner
# @affiliation: Medical Systems Biology Group, University of Lübeck
# @date: 2019-10-21
# @modified by Karin Näsvall 200326
###############

# R >= 3.6.x
# load libraries
#
#on Uppmax module load RStudio/1.1.463  r_packages/3.6.1
library(phyloseq)
library(Biostrings)
library(gdata)
library(dplyr)
library(ggplot2)
library(nlme)
library(ade4)
library(ape)
library(vegan)
library(gridExtra)
library(reshape2)
library(wesanderson)
library(viridis)
library(scales)

result_list <- list()

#
# Read data from disk
# data is located in a folder called 'data'
#
# 1. ASV table from dada2
# 2. taxonomic classification (from dada2/DECIPHER)
# 3. read phylogenetic tree
#

seqtab.nochim <- readRDS("../data/FW_final_seqtab.rds") # CHANGE ME ...
taxid <- readRDS("../data/FW_final_tax.rds") # CHANGE ME ...
fitGTR <- readRDS("../data/FW_final_tree.rds") # CHANGE ME ...

# 4. read metadata 

#changed the name of teh levels in pop and hostplant in meta2
meta <- read.delim("../data/meta2.txt")
meta$X.ID <- gsub("P12708_", "", meta$X.ID) # clip prefix
rownames(meta) <- meta$X.ID; meta$X.ID <- NULL



# create phyloseq object 'ps',
# rename ASVs and add sequence information to phyloseq object

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(meta), 
               tax_table(taxid),
               phy_tree(fitGTR$tree))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
sample_data(ps)$ID <- sample_names(ps)
ps

#
# Check taxa distribution and keep only bacteria and remove:
# Phylum not assigned, Mitochondria,Chloroplast
#

table( tax_table(ps)[, "Kingdom"], useNA = "always")
table( tax_table(ps)[, "Phylum"], useNA = "always")
table( tax_table(ps)[, "Class"], useNA = "always")
table( tax_table(ps)[, "Order"], useNA = "always") 
table( tax_table(ps)[, "Family"], useNA = "always") # remove Mitochondria
table( tax_table(ps)[, "Genus"], useNA = "always") # remove Chloroplast

# remove...
x <- subset_taxa(ps, Kingdom=="Bacteria")
x <- subset_taxa(x, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized", "unclassified"))
x <- subset_taxa(x, Genus!="Chloroplast")
x <- subset_taxa(x, Family!="Mitochondria")

x_wolbachia <- subset_taxa(x, Genus=="Wolbachia")
x_wolbachia

#
# check yield after preprocessing
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 42.0    693.8   2431.0  17684.3   9233.5 199575.0 
# > sd( samsum ); mad( samsum )
# [1] 36630.24
# [1] 3134.216
samsum <- sort( sample_sums(x) )
summary( samsum )
sd( samsum ); mad( samsum )

pdf(file = "plots/contigs.pdf", height = 6, width = 6, useDingbats = F)
par(las=1)
plot(log10(samsum), pch=19, col="skyblue3", 
     main = "Contigs per sample after preprocessing/filtering", ylab="#Contigs (*1000)", xlab="Index", axes=F, ylim=c(0,6))
axis(1, tick = F)
axis(2, at=c(0, 2, 4, 6), labels = c(1, 0.1, 10, 100))
sort( apply ( data.frame(otu_table(x)), 1, sum ) )
abline(h=log10(500), col="grey45", lty=2)
dev.off()

#
# remove samples < 500 contigs (prune samples)
# and root phylogenetic tree
#

x2 <- prune_samples(sample_sums(x) >= 500, x)
# outlier determined using weighted UniFrac and taxa distribution
outliers <- c("1024") 
x2 <- prune_samples(!(sample_names(x2) %in% outliers), x2)
# remove singletons, might be just sequencing errors
fil <- genefilter_sample( x2, filterfun_sample(function(x) x >= 1), A = 1 ); sum(fil)
( x2 <- prune_taxa(fil, x2) )

# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 178 taxa and 54 samples ]
# sample_data() Sample Data:       [ 54 samples by 6 sample variables ]
# tax_table()   Taxonomy Table:    [ 178 taxa by 7 taxonomic ranks ]
# phy_tree()    Phylogenetic Tree: [ 178 tips and 176 internal nodes ]
# refseq()      DNAStringSet:      [ 178 reference sequences ]

# set root
ape::is.rooted(phy_tree(x2))
phy_tree(x2) = ape::root(phy_tree(x2), sample(taxa_names(x2), 1), resolve.root = TRUE)
is.rooted(phy_tree(x2))

x2_wolbachia <- subset_taxa(x2, Genus=="Wolbachia")
x2_wolbachia

#
# rarefy to 500 contigs per sample
# and log transform abundance table
# for Unifrac plots (note: log transformation not necessary)
#

x.500 <- rarefy_even_depth(x2, sample.size = 500, replace = FALSE, rngseed = 123)
x.500.log <- transform_sample_counts(x.500, function(x) log(1 + x))

# 75OTUs were removed because they are no longer 
# present in any sample after random subsampling

#
# check Beta diversity
# and plot top 20 Genera
# 

b1 <- plot_ordination(x.500.log, 
                ordinate(x.500, method="PCoA", distance="bray"), 
                color="Population", shape = "Hostplant", label = "ID", 
                title="Bray PCoA")
b2 <- plot_ordination(x.500.log, 
                ordinate(x.500, method="PCoA", distance="unifrac"), 
                color="Population", shape = "Hostplant", label = "ID", 
                title="Unifrac PCoA")
b3 <- plot_ordination(x.500.log, 
                      ordinate(x.500, method="PCoA", distance="wunifrac"), 
                      color="Population", shape = "Hostplant", label = "ID", 
                      title="w. Unifrac PCoA")

gridExtra::grid.arrange(b2, b3, b1, ncol = 2)

top20 <- names(sort(taxa_sums(x2), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(x2, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 

#top20 without Wolbachia
x2_noW <- subset_taxa(x2, Genus != "Wolbachia")
top20_noW <- names(sort(taxa_sums(x2_noW), decreasing=TRUE))[1:20]
ps.top20_noW <- transform_sample_counts(x2, function(OTU) OTU/sum(OTU))
ps.top20_noW <- prune_taxa(top20_noW, ps.top20_noW)
plot_bar(ps.top20_noW, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 

plot_bar(x2, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 


top50 <- names(sort(taxa_sums(x2), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(x2, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)


#plot nr of reads per sample
pdf(file = "plots/nr_of_reads_x2.pdf", height = 6, width = 6, useDingbats = F)
plot(sort(sample_sums(x2), TRUE), type = "h", ylab = "reads" 
    )
dev.off()

pdf(file = "plots/nr_of_reads_x500.pdf", height = 6, width = 6, useDingbats = F)
plot(sort(sample_sums(x.500), TRUE), type = "h", ylab = "reads" 
)
dev.off()

#merge samples by population to make a plot of percentages
#OBS merging sums the data from the merged samples so do not look at absolute values

#first transform, ie normalising for different read depth in the different samples
x2_tr = transform_sample_counts(x2, function(a) 100 * a/sum(a))
#then merge to account for different numbers of samples in each group
x2_tr.merged = x2_tr

variable1 = as.character(get_variable(x2_tr.merged, "Population"))
variable2 = as.character(get_variable(x2_tr.merged, "Hostplant"))
sample_data(x2_tr.merged)$X_F <- mapply(paste0, variable1, variable2, 
                                              collapse = "_")
x2_tr.merged = merge_samples(x2_tr.merged, "X_F")

#transform again
x2_tr.merged.100  = transform_sample_counts(x2_tr.merged, function(x) 100 * x/sum(x))
plot_bar(x2_tr.merged.100, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x")


ps.top20_x2 <- prune_taxa(top20, x2_tr.merged.100)
plot_bar(ps.top20_x2, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 

#merge top 20
#first transform, ie normalising for different read depth in the different samples
#ps.top20_tr = transform_sample_counts(ps.top20, function(a) 100 * a/sum(a))
#then merge to account for different numbers of samples in each group
#ps.top20_tr.merged = ps.top20_tr

# 
# x2_merged <- merge_samples(x2, c("Population", "Hostplant"))
# sample_data(x2_merged)$Population <- levels(sample_data(x2)$Population)
# #get percentages
# x2_merged = transform_sample_counts(x2_merged, function(a) 100 * a/sum(a))
# #prune all but the top20 taxa
# x2_merged_top20 <- prune_taxa(top20, x2_merged)
# 
# 

#change order default is alphabetical, get the same order an colour as genus levels
# p$data$Order <- as.character(p$data$Order)
# p$data$Order <- factor(x = p$data$Order, 
#                        levels = unique(p$data$Order[order(as.character(p$data$Phylum))]))


##############################################################################

##############################################################################

#Tax level
#cat
#top25
top25cat <- names(sort(taxa_sums(subset_samples(x2, sample_data(x2)$Population == "Catalonia")), decreasing=TRUE))[1:25]
ps.top25_cat <- transform_sample_counts(x2, function(OTU) OTU/sum(OTU))
#using these top 25 on the merged dataset 
ps.top25_cat <- prune_taxa(top25cat, x2_tr.merged.100)
ps.top25_cat <- subset_samples(ps.top25_cat, sample_data(ps.top25_cat)$Population == "1")
#do not need to tranform again could be good to see how much has been pruned away?
#ps.top25_cat = transform_sample_counts(ps.top25_cat, function(a) 100 * a/sum(a))
plot_bar(ps.top25_cat, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 

#change order of level hostplant is in alphabetical order... 1 lotus 2 dorycnium
#p$data$Order <- factor(x = p$data$Order, 
#                       levels =
sample_data(ps.top25_cat)$Hostplant <- factor(x=sample_data(ps.top25_cat)$Hostplant, levels = c("2","1"))
sample_data(ps.top25_cat)$Hostplant
#[1] 1 2
#Levels: 2 1

#for phylum
sample_data(x2_tr.merged.100)$Hostplant <- factor(x=sample_data(x2_tr.merged.100)$Hostplant, levels = c("2","1"))
sample_data(x2_tr.merged.100)$Hostplant


#in one figure with grid.arrange

#colour spec
show_col(viridis_pal()(72))
show_col(magma_pal()(72))


#MyColour <- c("#FF0000", "#00FF00", "#0000FF") 
#names(MyColour) <- c("Mazda RX4", "Toyota Corolla", "Fiat 128")

col_cat_genus <- c("Wolbachia" = "#287D8EFF", "Kocuria" = "#55C667FF", "Burkholderia" = "#39568CFF",  "Kineococcus"= "#404788FF",   "Methylobacterium" ="#95D840FF", "Pseudomonas_E" = "#B8DE29FF",  "Rothia" = "#DCE319FF", "Sphingomonas" = "#FDE725FF")         

col_cat_family<- c("#287D8EFF", "#55C667FF", "#404788FF", "#B8DE29FF", "#482677FF", "#39568CFF", "#FDE725FF")
names(col_cat_family)<-get_taxa_unique(ps.top25_cat, "Family")

col_cat_phylum_viridis<-
    c("#287D8EFF",
      "#55C667FF",
      "#440154FF",
      "#481567FF",
      "#482677FF",
      "#453781FF",
      "#404788FF",
      "#39568CFF",
      "#95D840FF",
      "#B8DE29FF",
      "#DCE319FF",
      "#E3E418FF",
      "#FDE725FF")

#magma
col_cat_phylum<-
    c("#287D8EFF",
      "#55C667FF",
      "#FCFFB2",
      "#FBC17D",
      "#FA8657",
      "#ED504A",
      "#E03B50",
      "#C92D59",
      "#B02363",
      "#981D69",
      "#81176D",
      "#6B116F",
      "#57096E")


get_taxa_unique(x2_tr.merged.100, "Phylum")
#sorted in order
names(col_cat_phylum)<-
    c("Proteobacteria",    
      "Actinobacteriota", 
      "Acidobacteriota", 
      "Bacteroidota",   
      "Campylobacterota", 
      "Fibrobacterota",  
      "Firmicutes",   
      "Firmicutes_A", 
      "Firmicutes_C",  
      "Firmicutes_I",
      "Fusobacteriota",
      "Myxococcota",
      "Verrucomicrobiota")

    
plot_cat_genus<-plot_bar(ps.top25_cat, fill = "Genus", x = "Hostplant", y = "Abundance", title = " ") +
    ylab("Proportion of Sequences") +
    geom_bar(aes(fill=Genus), stat = "identity", position = "stack") + 
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, vjust = 1),
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(name=" ", labels=c("L. dorycnium", "L. corniculatus")) +
    scale_fill_manual(values = col_cat_genus)

#change order of levels for family
plot_family_cat <- plot_bar(ps.top25_cat, fill = "Family", x = "Hostplant", y = "Abundance", title = " ")
plot_family_cat$data$Family <- factor(x=plot_family_cat$data$Family, levels = unique(plot_family_cat$data$Family[order(as.character(plot_family_cat$data$Genus))]))
plot_cat_family<-plot_family_cat +
    ylab("Proportion of Sequences") +
    geom_bar(aes(fill=Family), stat = "identity", position = "stack") + 
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, vjust = 1),
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(name=" ", labels=c("L. dorycnium", "L. corniculatus")) +
    scale_fill_manual(values = col_cat_family)


#plot Phylum with the whole dataset
title="Catalonia"
plot_cat_phylum<-plot_bar(subset_samples(x2_tr.merged.100, sample_data(x2_tr.merged.100)$Population == "1"), fill = "Phylum", x = "Hostplant", y = "Abundance", title = title) +
    ylab("Proportion of Sequences") +
    geom_bar(aes(fill=Phylum), stat = "identity", position = "stack") + 
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, vjust = 1),
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_x_discrete(name=" ", labels=c("L. dorycnium", "L. corniculatus")) +
    scale_fill_manual(values = col_cat_phylum)    

grid.arrange(plot_cat_phylum, plot_cat_family, plot_cat_genus, nrow=1)

##############
#SWE
#top25
top25swe <- names(sort(taxa_sums(subset_samples(x2, sample_data(x2)$Population == "Sweden")), decreasing=TRUE))[1:25]
ps.top25_swe <- transform_sample_counts(x2, function(OTU) OTU/sum(OTU))
#using these top 25 on the merged dataset 
ps.top25_swe <- prune_taxa(top25swe, x2_tr.merged.100)
ps.top25_swe <- subset_samples(ps.top25_swe, sample_data(ps.top25_swe)$Population == "2")
#do not need to tranform again could be good to see how much has been pruned away?
#ps.top25_swe = transform_sample_counts(ps.top25_swe, function(a) 100 * a/sum(a))
plot_bar(ps.top25_swe, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 

sample_data(ps.top25_swe)$Hostplant <- factor(x=sample_data(ps.top25_swe)$Hostplant, levels = c("2","1"))
sample_data(ps.top25_swe)$Hostplant


#in one figure with grid.arrange
#order colour
#get_taxa_unique(ps.top25_cat, "Family")
col_swe_genus <- c("#287D8EFF", "#55C667FF", "#39568CFF", "#453781FF", "#404788FF", "#95D840FF", "#B8DE29FF", "#DCE319FF", "#E3E418FF", "#FDE725FF")
names(col_swe_genus)<- c(
    "Wolbachia", 
    "Kocuria", 
    "Burkholderia",
    "Cloacibacterium",
    "Cutibacterium",
    "Methylobacterium",
    "Psychrobacter",   
    "Sphingomonas",
    "Staphylococcus",
    "Stenotrophomonas")


plot_swe_genus<-plot_bar(ps.top25_swe, fill = "Genus", x = "Hostplant", y = "Abundance", title = " ") +
    ylab("Proportion of Sequences") +
    geom_bar(aes(fill=Genus), stat = "identity", position = "stack") + 
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, vjust = 1),
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(name=" ", labels=c("L. dorycnium", "L. corniculatus")) +
    scale_fill_manual(values = col_swe_genus)
#order colours
col_swe_family <- c("#287D8EFF", 
                    "#55C667FF", 
                    "#39568CFF", 
                    "#453781FF", 
                    "#404788FF", 
                    "#95D840FF", 
                    "#B8DE29FF",
                    "#DCE319FF",
                    "#E3E418FF",
                    "#FDE725FF")


names(col_swe_family)<- c("Anaplasmataceae", 
                          "Micrococcaceae",
                          "Burkholderiaceae",
                          "Beijerinckiaceae",
                          "Moraxellaceae",
                          "Propionibacteriaceae",     
                          "Xanthomonadaceae",               
                          "Staphylococcaceae",
                          "Weeksellaceae",
                          "Sphingomonadaceae")


#change order of levels for family
plot_family_swe <- plot_bar(ps.top25_swe, fill = "Family", x = "Hostplant", y = "Abundance", title = " ")
plot_family_swe$data$Family <- factor(x=plot_family_swe$data$Family, levels = unique(plot_family_swe$data$Family[order(as.character(plot_family_swe$data$Genus))]))
plot_swe_family<-plot_family_swe +
    ylab("Proportion of Sequences") +
    geom_bar(aes(fill=Family), stat = "identity", position = "stack") + 
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, vjust = 1),
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    scale_x_discrete(name=" ", labels=c("L. dorycnium", "L. corniculatus")) +
    scale_fill_manual(values = col_swe_family)

#change colour order, the number of taxa
#in cat and swe are the same on phylum level
#magma
col_swe_phylum<-
    c("#287D8EFF",
      "#55C667FF",
      "#FCFFB2",
      "#FBC17D",
      "#FA8657",
      "#ED504A",
      "#E03B50",
      "#C92D59",
      "#B02363",
      "#981D69",
      "#81176D",
      "#6B116F",
      "#57096E")


get_taxa_unique(subset_samples(x2, sample_data(x2)$Population == "Sweden"), "Phylum")
#sorted in order
names(col_swe_phylum)<-
    c("Proteobacteria",    
      "Actinobacteriota", 
      "Acidobacteriota", 
      "Bacteroidota",   
      "Campylobacterota", 
      "Fibrobacterota",  
      "Firmicutes",   
      "Firmicutes_A", 
      "Firmicutes_C",  
      "Firmicutes_I",
      "Fusobacteriota",
      "Myxococcota",
      "Verrucomicrobiota")

col_swe_phylum_top10<-
    c("#287D8EFF",
      "#55C667FF",
      "#FCFFB2",
      "#FBC17D",
      "#FA8657",
      "#ED504A",
      "#E03B50",
      "#C92D59",
      "#B02363",
      "#981D69")
names(col_swe_phylum_top10)<- c("Proteobacteria",    
                                   "Actinobacteriota", 
                                   "Acidobacteriota", 
                                   "Bacteroidota", 
                                   "Fibrobacterota",  
                                   "Firmicutes",   
                                   "Firmicutes_A", 
                                   "Firmicutes_C",  
                                   "Firmicutes_I",
                                   "Fusobacteriota")

    

#plot Phylum with the whole dataset
title="Sweden"
plot_swe_phylum<-plot_bar(subset_samples(x2_tr.merged.100, sample_data(x2_tr.merged.100)$Population == "2"), fill = "Phylum", x = "Hostplant", y = "Abundance", title = title) +
    ylab("Proportion of Sequences") +
    geom_bar(aes(fill=Phylum), stat = "identity", position = "stack") + 
    theme(panel.border = element_blank(), 
          axis.line = element_line(color = 'black'),
          axis.text.x = element_text(angle = 30, face = "italic", hjust = 1, vjust = 1),
          #axis.text.y = element_text(angle = 90, hjust = 0.5),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())+
    scale_x_discrete(name=" ", labels=c("L. dorycnium", "L. corniculatus")) +
    scale_fill_manual(values = col_swe_phylum)    


pdf(file = "plots/tax_lev.pdf", useDingbats = F)


grid.arrange(plot_cat_phylum, plot_cat_family, plot_cat_genus, plot_swe_phylum, plot_swe_family, plot_swe_genus, nrow=2)


#get tables with relative abundance
rel_abundance_phylum.table <- otu_table(tax_glom(x2_tr.merged.100, "Phylum"))
rel_abundance_phylum.table <- as.data.frame(t(rel_abundance_phylum.table))

rel_abundance_family.table <- otu_table(tax_glom(x2_tr.merged.100, "Family"))
rel_abundance_family.table <- as.data.frame(t(rel_abundance_family.table))

rel_abundance.table <- otu_table(tax_glom(x2_tr.merged.100, "Genus"))
rel_abundance.table <- as.data.frame(t(rel_abundance.table))

#get taxa, genus level
tax_table(x2_tr.merged.100)[,6]

#test
#Phylum
#cat
wilcox.test(rel_abundance_phylum.table[,1], rel_abundance_phylum.table[,2])
#swe
wilcox.test(rel_abundance_phylum.table[,3], rel_abundance_phylum.table[,4])

#Family
#cat
wilcox.test(rel_abundance_family.table[,1], rel_abundance_family.table[,2])
#swe
wilcox.test(rel_abundance_family.table[,3], rel_abundance_family.table[,4])

#Genus
#cat
wilcox.test(rel_abundance.table[,1], rel_abundance.table[,2])
#swe
wilcox.test(rel_abundance.table[,3], rel_abundance.table[,4])

#testing taxa level
#mt(physeq, classlabel, minPmaxT = "minP", method = "fdr", ...)
#This is a function for multiple testing of difference in taxonomic abundance with Wilcoxon ranked test 
#betweeen variable in phyloseq objects, I run with method = "fdr", default Benjamini-Hoch correction for multiple testing

#test
result_mt_pop <- mt(x2_tr, "Population", method = "BH")
result_mt_host <- mt(x2_tr, "Hostplant", method = "BH")

#each pop
subset_samples(x2_tr, sample_data(x2_tr)$Population == "Catalonia")
subset_samples(x2_tr, sample_data(x2_tr)$Population == "Sweden")

result_mt_cat_genus <- mt(subset_samples(tax_glom(x2_tr, "Genus"), 
                                         sample_data(tax_glom(x2_tr, "Genus"))$Population == "Catalonia"), 
                          "Hostplant", 
                          method = "fdr", 
                          test = "wilcoxon", 
                          nonpara = "y")

result_mt_swe_genus <- mt(subset_samples(tax_glom(x2_tr, "Genus"), 
                                         sample_data(tax_glom(x2_tr, "Genus"))$Population == "Sweden"), 
                          "Hostplant", 
                          method = "fdr", 
                          test = "wilcoxon", 
                          nonpara = "y")

#family
result_mt_cat_family <- mt(subset_samples(tax_glom(x2_tr, "Family"), 
                                         sample_data(tax_glom(x2_tr, "Family"))$Population == "Catalonia"), 
                          "Hostplant", 
                          method = "fdr", 
                          test = "wilcoxon", 
                          nonpara = "y")

result_mt_swe_family <- mt(subset_samples(tax_glom(x2_tr, "Family"), 
                                         sample_data(tax_glom(x2_tr, "Family"))$Population == "Sweden"), 
                          "Hostplant", 
                          method = "fdr", 
                          test = "wilcoxon", 
                          nonpara = "y")

#phylum
result_mt_cat_phylum <- mt(subset_samples(tax_glom(x2_tr, "Phylum"), 
                                          sample_data(tax_glom(x2_tr, "Phylum"))$Population == "Catalonia"), 
                           "Hostplant", 
                           method = "fdr", 
                           test = "wilcoxon", 
                           nonpara = "y")

result_mt_swe_phylum <- mt(subset_samples(tax_glom(x2_tr, "Phylum"), 
                                          sample_data(tax_glom(x2_tr, "Phylum"))$Population == "Sweden"), 
                           "Hostplant", 
                           method = "fdr", 
                           test = "wilcoxon", 
                           nonpara = "y")




write.table(rel_abundance.table, file = "rel_abundance_genus.table")
write.table(rel_abundance_family.table, file = "rel_abundance_family.table")
write.table(rel_abundance_phylum.table, file = "rel_abundance_phylum.table")

write.table(result_mt_cat_genus, file = "diff_tax_abundance_mt_test_genus_cat.table")
write.table(result_mt_swe_genus, file = "diff_tax_abundance_mt_test_genus_swe.table")
write.table(result_mt_cat_family, file = "diff_tax_abundance_mt_test_family_cat.table")
write.table(result_mt_swe_family, file = "diff_tax_abundance_mt_test_family_swe.table")
write.table(result_mt_cat_phylum, file = "diff_tax_abundance_mt_test_phylum_cat.table")
write.table(result_mt_swe_phylum, file = "diff_tax_abundance_mt_test_phylum_swe.table")



#####################
#without Wolbachia
#cat
#for labelling facets in ggplot
hostplant_names <- c("2"="L.dorycnium", "1"="L.corniculatus")

x2_noW <- subset_taxa(x2, Genus != "Wolbachia")
x2_tr.merged_noW <- subset_taxa(x2_tr.merged, Genus != "Wolbachia")
#x2_tr.merged.100_noW  = transform_sample_counts(x2_tr.merged_noW, function(x) 100 * x/sum(x))

top25cat_noW <- names(sort(taxa_sums(subset_samples(x2_noW, sample_data(x2_noW)$Population == "Catalonia")), decreasing=TRUE))[1:25]
ps.top25_cat_noW <- transform_sample_counts(x2_noW, function(OTU) OTU/sum(OTU))
#using these top 25 on the merged dataset 
ps.top25_cat_noW <- prune_taxa(top25cat_noW, x2_tr.merged_noW)
#ps.top25_cat_noW <- subset_samples(ps.top25_cat_noW, sample_data(ps.top25_cat_noW)$Population == "1")
#do not need to tranform again could be good to see how much has been pruned away?
#ps.top25_cat = transform_sample_counts(ps.top25_cat, function(a) 100 * a/sum(a))
#plot_bar(ps.top25_cat, fill="Genus") + facet_wrap(~Population+Hostplant, scales="free_x") 

#change order of level hostplant is in alphabetical order... 1 lotus 2 dorycnium
#p$data$Order <- factor(x = p$data$Order, 
#                       levels =
sample_data(ps.top25_cat_noW)$Hostplant <- factor(x=sample_data(ps.top25_cat_noW)$Hostplant, levels = c("2","1"))
sample_data(ps.top25_cat_noW)$Hostplant
#[1] 1 2
#Levels: 2 1

top25swe_noW <- names(sort(taxa_sums(subset_samples(x2_noW, sample_data(x2_noW)$Population == "Sweden")), decreasing=TRUE))[1:25]
ps.top25_swe_noW <- transform_sample_counts(x2_noW, function(OTU) OTU/sum(OTU))
#using these top 25 on the merged dataset 
ps.top25_swe_noW <- prune_taxa(top25swe_noW, x2_tr.merged_noW)
sample_data(ps.top25_swe_noW)$Hostplant <- factor(x=sample_data(ps.top25_swe_noW)$Hostplant, levels = c("2","1"))
sample_data(ps.top25_swe_noW)$Hostplant

col_cat_genus_woW <- c("#55C667FF",
                       "#440154FF",
                       "#481567FF",
                       "#482677FF",
                       "#453781FF",
                       "#404788FF",
                       "#39568CFF",
                       "#33638DFF",
                       "#2D708EFF",
                       "#238A8DFF",
                       "#20A387FF",
                       "#B8DE29FF",
                       "#DCE319FF",
                       "#FDE725FF")


names(col_cat_genus_woW)<- c("Kocuria", 
  "Burkholderia",
  "Corynebacterium", 
  "Cutibacterium",
  "Kineococcus",
  "Leaf454",
  "Methylobacterium",        
  "Micrococcus",           
  "Moraxella_A",
  "Nocardioides",                  
  "Pseudomonas_E",                  
  "Rothia",
  "Sphingomonas",
  "Stenotrophomonas")

col_swe_genus_woW <- c("#55C667FF",
                       "#460B5EFF",
                       "#481668FF",
                       "#482070FF",
                       "#440154FF",
                       "#472D7BFF",
                       "#481567FF",
                       "#482677FF",
                       "#453781FF",
                       "#2D708EFF",
                       "#404788FF",
                       "#39568CFF",
                       "#33638DFF",
                       "#D0E11CFF",
                       "#DCE319FF",
                       "#F5E61FFF",
                       "#FDE725FF")

names(col_swe_genus_woW) <- c("Kocuria",
    "Acinetobacter",
    "Anoxybacillus",
    "Bacillus_E",
    "Burkholderia",
    "Cloacibacterium",  
    "Corynebacterium",
    "Cutibacterium",
    "Herbaspirillum",
    "Lawsonella",
    "Leaf454",          
    "Methylobacterium",
    "Micrococcus",
    "Psychrobacter",        
    "Sphingomonas",
    "Staphylococcus",
    "Stenotrophomonas")

pdf(file = "plots/tax_lev_cat_noWtop25.pdf", height = 6, width = 6, useDingbats = F)
plot_cat_woWol<-plot_bar(ps.top25_cat_noW, "Genus", fill = "Genus", title = "Catalonia") + 
    facet_grid(.~Hostplant, labeller = labeller(Hostplant = hostplant_names)) +
    geom_bar(aes(fill=Genus), stat = "identity", position = "stack") +
    ylab("Proportion of Sequences") +
    xlab("") + 
    ylim(0,30) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(face = "italic", size = 12)) +
    scale_fill_manual(values=col_cat_genus_woW, 
                      guide=
                          guide_legend(ncol = 2))    

dev.off()

pdf(file = "plots/tax_lev_swe_noWtop25.pdf", height = 6, width = 6, useDingbats = F)
plot_swe_woWol<-plot_bar(ps.top25_swe_noW, "Genus", fill = "Genus", title = "Sweden") + 
    facet_grid(.~Hostplant, labeller = labeller(Hostplant = hostplant_names)) +
    geom_bar(aes(fill=Genus), stat = "identity", position = "stack") +
    ylab("Proportion of Sequences") +
    xlab("") + 
    ylim(0,30) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.ticks = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.text = element_text(face = "italic", size = 12)) +
    scale_fill_manual(values = viridis(17), 
                      guide=
                          guide_legend(ncol = 2))    

dev.off()

grid.arrange(plot_cat_woWol, plot_swe_woWol)


##########################################################################

###########################################################################

#testing a tree
pdf(file = "plots/tree_test.pdf", height = 6, width = 6, useDingbats = F)
plot_tree(ps.top20, color=c("Population"), label.tips="Genus", ladderize="left", plot.margin=0.3)
dev.off()

pdf(file = "plots/tree_test_hostplant.pdf", height = 6, width = 6, useDingbats = F)
plot_tree(ps.top20, color="Hostplant", label.tips="Genus", ladderize="left", plot.margin=0.3)
dev.off()

pdf(file = "plots/tree_test_ASV.pdf", height = 6, width = 6, useDingbats = F)
plot_tree(x_wolbachia, color=c("Population"), label.tips="taxa_names", ladderize="left", plot.margin=0.3)
dev.off()


##################################################
# Alpha diversity; Shannon index
##################################################


ps_alpha_div <- estimate_richness(x2, split = TRUE, measure = "Shannon")
ps_alpha_div$ID <- rownames(ps_alpha_div) %>%
    as.factor()
ps_samp <- sample_data(x2) %>%
    unclass() %>%
    data.frame() %>%
    left_join(ps_alpha_div, by = "ID") %>%
    reshape2::melt(measure.vars = "Shannon",
         variable.name = "diversity_measure",
         value.name = "alpha_diversity")

# reorder's facet from lowest to highest diversity
diversity_means <- ps_samp %>%
    group_by(Population) %>%
    summarise(mean_div = mean(alpha_diversity)) %>%
    arrange(mean_div)
ps_samp$Hostplant <- factor(ps_samp$Hostplant)

# Linear mixed model for Population
alpha_div_model_pop <- lme(fixed = alpha_diversity ~ Population, data = ps_samp,
                           random = ~ 1 | Hostplant)

new_data <- expand.grid(Hostplant = levels(ps_samp$Hostplant),
                        Population = levels(ps_samp$Population))
new_data$pred <- predict(alpha_div_model_pop, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model_pop$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model_pop$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model_pop$sigma ^ 2

# fitted values, with error bars
pdf(file = "plots/Alpha_population.pdf", height = 6, width = 6, useDingbats = F)
ggplot(ps_samp %>% left_join(new_data)) +
    geom_errorbar(aes(x = Population, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(x = Population, y = alpha_diversity,
                   col = Population), size = 2.5) +
    facet_wrap(~Hostplant) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    scale_color_manual(values = c("deepskyblue2", "orange1")) +
    
    labs(x = "", y = "Shannon Diversity", color = "Population") +
    guides(col = guide_legend(override.aes = list(size = 5))) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
dev.off()

#new colours for manuscript, filled circles
ps_samp_new <- ps_samp %>% left_join(new_data)
levels(ps_samp_new$Hostplant) <- factor(x=ps_samp_new$Hostplant, levels = c("L. dorycnium", "L. corniculatus"))

ggplot(ps_samp_new, aes(Population)) +
    geom_errorbar(aes(Population, alpha_diversity, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(Population, alpha_diversity, 
                   color=Population, 
                   shape=Hostplant, 
                   fill=interaction(Population,Hostplant)), size = 3, stroke = 1) +
    facet_wrap(~Hostplant) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    scale_shape_manual(values = c(19, 21),
                       guide = guide_legend(label.theme = element_text(angle = 0, face = "italic"))) +
    scale_color_manual(values = c("orangered4", "orange3")) +
    scale_fill_manual(values = c("orangered4", "orange1", "tomato1", "gold"),
                      guide = FALSE) +
    labs(x = "", y = "Shannon Diversity") +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks.x = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 12),
          strip.text = element_text(face = "italic", size=18))




# Linear mixed model for Hostplant
alpha_div_model_host <- lme(fixed = alpha_diversity ~ Hostplant, data = ps_samp,
                            random = ~ 1 | Population)

new_data <- expand.grid(Population = levels(ps_samp$Population),
                        Hostplant = levels(ps_samp$Hostplant))
new_data$pred <- predict(alpha_div_model_host, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model_host$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model_host$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model_host$sigma ^ 2


# fitted values, with error bars
pdf(file = "plots/Alpha_hostplant.pdf", height = 6, width = 6, useDingbats = F)
ggplot(ps_samp %>% left_join(new_data)) +
    geom_errorbar(aes(x = Hostplant, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(x = Hostplant, y = alpha_diversity,
                   col = Hostplant), size = 2.5) +
    facet_wrap(~Population) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    scale_color_brewer(palette = "Set2") +
    labs(x = "", y = "Shannon Diversity", color = "Hostplant") +
    guides(col = guide_legend(override.aes = list(size = 5))) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
dev.off()

ps_samp_new <- ps_samp %>% left_join(new_data)
#change order of levels
ps_samp_new$Hostplant <- factor(x=ps_samp_new$Hostplant, levels = c("L.dorycnium", "L.corniculatus"))
#change name
levels(ps_samp_new$Hostplant)<- c("L. dorycnium", "L. corniculatus")

#green
#ggplot(ps_samp_new, aes(Hostplant)) +
    geom_errorbar(aes(Hostplant, alpha_diversity, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(Hostplant, alpha_diversity, 
                   color=Hostplant, 
                   shape=Population), 
                   #fill=interaction(Population,Hostplant)),
                   size = 3, stroke = 2) +
    facet_wrap(~Population) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    scale_shape_manual(values = c(19, 21),
                       guide = FALSE) +
    scale_color_manual(values = c("chartreuse4", "chartreuse2"), 
                       guide = guide_legend(label.theme = element_text(angle = 0,
                                                                       face = "italic", 
                                                                       size = 16))) +
    #scale_fill_manual(values = c("chartreuse4","chartreuse4", "chartreuse2","chartreuse2"),
    #                  guide = FALSE) +
    labs(x = "", y = "Shannon Diversity") +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0), size = 2),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks.x = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          strip.text = element_text(size=20), 
          strip.background = element_rect(color = "black", fill = "white", size = 2))

#with filled dots no pop diff
ggplot(ps_samp_new, aes(Hostplant)) +
    geom_errorbar(aes(Hostplant, alpha_diversity, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(Hostplant, alpha_diversity, 
                   fill=Hostplant), 
                   shape=21, 
               #fill=interaction(Population,Hostplant)),
               size = 3) +
    facet_wrap(~Population) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    #scale_shape_manual(values = c(19, 21),
    #                   guide = FALSE) +
    scale_fill_manual(values = c("chartreuse4", "chartreuse2"), 
                       guide = guide_legend(override.aes = list(size = 5, shape = 21),
                                            label.theme = element_text(angle = 0,
                                                                       face = "italic", 
                                                                       size = 16))) +
    #scale_fill_manual(values = c("chartreuse4","chartreuse4", "chartreuse2","chartreuse2"),
    #                  guide = FALSE) +
    labs(x = "", y = "Shannon Diversity") +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0), size = 2),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.ticks.x = element_blank(),
          legend.key = element_blank(),
          legend.text = element_text(size = 16),
          legend.title = element_text(size = 20),
          strip.text = element_text(size=20), 
          strip.background = element_rect(color = "black", fill = "white", size = 2))

# 
# Linear mixed model Instar
alpha_div_model_instar <- lme(fixed = alpha_diversity ~ Instar, 
                              data = ps_samp,
                              random = list(~ 1|Population, ~ 1|Hostplant) )

new_data <- expand.grid(Population = levels(ps_samp$Population),
                        Hostplant = levels(ps_samp$Hostplant),
                        Instar = levels(ps_samp$Instar))
new_data$pred <- predict(alpha_div_model_instar, newdata = new_data)
X <- model.matrix(eval(eval(alpha_div_model_instar$call$fixed)[-2]),
                  new_data[-ncol(new_data)])
pred_var_fixed <- diag(X %*% alpha_div_model_instar$varFix %*% t(X))
new_data$pred_var <- pred_var_fixed + alpha_div_model_instar$sigma ^ 2

# fitted values, with error bars
pdf(file = "plots/Alpha_Instar.pdf", height = 6, width = 6, useDingbats = F)
ggplot(ps_samp %>% left_join(new_data)) +
    geom_errorbar(aes(x = Instar, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(x = Instar, y = alpha_diversity,
                   col = Instar), size = 2.5) +
    facet_wrap(~Population+Hostplant) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    scale_color_manual(values = c("cornsilk4", "darkolivegreen3")) +
    
    labs(x = "", y = "Shannon Diversity", color = "Instar") +
    guides(col = guide_legend(override.aes = list(size = 5))) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
dev.off()


ps_samp_new2 <- ps_samp %>% left_join(new_data)
ps_samp_new2

ps_samp_new2$Hostplant <- factor(x=ps_samp_new2$Hostplant, levels = c("L.dorycnium", "L.corniculatus"))
levels(ps_samp_new2$Hostplant) <- c("L. dorycnium", "L. corniculatus")


#green with black contour
ggplot(ps_samp_new2) +
    geom_errorbar(aes(x = Instar, ymin = pred - 2 * sqrt(pred_var),
                      ymax = pred + 2 * sqrt(pred_var)),
                  col = "grey45", size = .5) +
    geom_point(aes(x = Instar, y = alpha_diversity,
                   shape = Instar, 
                   fill = Hostplant
                   ), 
               size = 4#, 
               #position = position_jitter(width=0.2, height = 0)
    ) +
    facet_grid(Hostplant~Population) +
    scale_y_continuous(limits = c(1, 3), breaks = seq(0, 5, .5)) +
    scale_fill_manual(values = c("chartreuse4", "chartreuse2")) +
    scale_shape_manual(values = c(21,25)) + 
    labs(x = "", y = "Shannon Diversity") +
    guides(fill = guide_legend(title = "Hostplant",
                              label.theme = element_text(size = 16, angle = 0, face = "italic"),
                              override.aes = list(size = 5, shape = 21)),
           shape = guide_legend(override.aes = list(size = 3, stroke = 2))) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0), size = 2),
          axis.text.x = element_text(angle = -90, size = 0),
          axis.text.y = element_text(size = 16),
          axis.title.y = element_text(size = 16),
          axis.ticks = element_blank(),
          strip.text.y = element_text(face="italic", size = 20),
          strip.text.x = element_text(size = 20),
          strip.background = element_rect(color = "black", fill = "white", size = 2),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key = element_rect(colour = NA, fill = NA),
          legend.text = element_text(size=16), 
          legend.title = element_text(size = 20)) 


#test boxplot
ggplot(ps_samp_new2, aes(x=Instar, y=alpha_diversity)) +
    geom_boxplot() +
    
#test boxplot
ggplot(ps_samp_new2, aes(x=Instar, y=alpha_diversity)) +
    geom_boxplot() +
    facet_grid(Population~Hostplant)

#test points with error bar alpha diversity
ggplot(ps_samp_new2, aes(x=Instar,y=alpha_diversity)) +
    facet_grid(Population~Hostplant) +
    geom_errorbar(aes(ymin=alpha_diversity + 1,  
                      ymax=alpha_diversity - 1))
    geom_point()

               

#summary(alpha_div_model_pop)
anova(alpha_div_model_pop)
#summary(alpha_div_model_host)
anova(alpha_div_model_host)
#summary(alpha_div_model_instar)
anova(alpha_div_model_instar)

# Traditional tests for population and hostplant: Alpha diversity Wilcoxon
wilcox.test( 
    ps_samp$alpha_diversity[ps_samp$Population == "Catalonia"],
    ps_samp$alpha_diversity[ps_samp$Population == "Sweden"])

wilcox.test( 
    ps_samp$alpha_diversity[ps_samp$Population == "Catalonia" & ps_samp$Hostplant == "L.dorycnium"],
    ps_samp$alpha_diversity[ps_samp$Population == "Catalonia" & ps_samp$Hostplant == "L.corniculatus"])

wilcox.test( 
    ps_samp$alpha_diversity[ps_samp$Population == "Sweden" & ps_samp$Hostplant == "L.dorycnium"],
    ps_samp$alpha_diversity[ps_samp$Population == "Sweden" & ps_samp$Hostplant == "L.corniculatus"])


#
# Beta diversity, here unweighted Unifrac due to structure of data
# (low number of contigs)
#


pdf(file = "plots/Beta_500.pdf", height = 6, width = 6, useDingbats = F)
plot_ordination(x.500.log, 
                ordinate(x.500.log, method="PCoA", "unifrac"), 
                type="sample", 
                shape="Hostplant", color="Population", #label="ID",
                title="", axes=c(1,2,3), justDF=F) + 
    scale_color_manual(values = c("deepskyblue2", "orange1")) +
    geom_point(size=4) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
dev.off()


pdf(file = "plots/Beta_500_manuscript.pdf", height = 6, width = 6, useDingbats = F)
plot_ordination(x.500.log, 
                ordinate(x.500.log, method="PCoA", "unifrac"), 
                type="sample",
                color = "Population",
                #label="ID",
                title="", axes=c(1,2,3), justDF=F) + 
    geom_point(aes(shape=Hostplant, color=Population, fill=Hostplant), stroke = 1.5, size=4) +    
   scale_shape_manual(values = c(19,21), 
                      guide = 
                            guide_legend(title = "Hostplant", label.theme = element_text(angle = 0, face = "italic", size = 10))) +
   scale_color_manual(values = c("orangered3", "orange1")) + 
   scale_fill_manual(values = c("transparent", "white"), 
                     guide = FALSE) +
   theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
         axis.text.x = element_text(angle = 0, size = 12),
         axis.text.y = element_text(size = 12),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(),
         legend.key = element_blank())
dev.off()

# PERMANOVA test
otu.uf  <- UniFrac(x.500.log, weighted=FALSE)
adonis(otu.uf ~ Population * Hostplant + Instar, 
       data=data.frame(sample_data(x.500.log)), 
       permutations = how(nperm=99999), parallel=8)
# Goodness-of-fit test
envfit( cmdscale(otu.uf) ~  Population, data=data.frame(sample_data(x.500.log)), permutations = 99999 )

# test beta for swedish and spanish (PERMANOVA)

x_spa <- subset_samples(x.500.log, Population == "Catalonia")
x_swe <- subset_samples(x.500.log, Population == "Sweden")

otu.uf  <- UniFrac(x_spa, weighted=FALSE)
adonis(otu.uf ~ Hostplant, 
       data=data.frame(sample_data(x_spa)), 
       permutations = how(nperm=99999), parallel=8)

otu.uf  <- UniFrac(x_swe, weighted=FALSE)
adonis(otu.uf ~ Hostplant, 
       data=data.frame(sample_data(x_swe)), 
       permutations = how(nperm=99999), parallel=8)

pdf(file = "plots/Beta_spanish.pdf", height = 6, width = 6, useDingbats = F)
plot_ordination(x_spa, 
                ordinate(x_spa, method="PCoA", "unifrac"), 
                type="sample", 
                color="Hostplant", #label="ID",
                title="", axes=c(1,2), justDF=F) + 
    #scale_color_manual(values = c("deepskyblue2", "orange1")) +
    scale_color_brewer(palette = "Set2") +
    geom_point(size=4) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank()) 
dev.off()
#change order of levels
sample_data(x_spa)$Hostplant <- factor(x=sample_data(x_spa)$Hostplant, levels = c("L.dorycnium","L.corniculatus"))
sample_data(x_spa)$Hostplant
levels(sample_data(x_spa)$Hostplant) <- c("L. dorycnium", "L. corniculatus")

#green
p <-plot_ordination(x_spa, 
                ordinate(x_spa, method="PCoA", "unifrac"), 
                type="sample", 
                #color="Hostplant",
                 #label="ID",
                title="Catalonia", axes=c(1,2), justDF=F#, 
                #position=position_jitter(seed=10, width = 0.5, height = 0)
) +
    geom_point(aes(fill=Hostplant),
               size=4,
               shape=21, 
               position=position_jitter(width = 0.05, height = 0.05, seed = 10)
               ) +
    #scale_shape_manual(values = c(21,19),
    scale_fill_manual(values = c("chartreuse4", "chartreuse2"),
                       guide=
                           guide_legend(override.aes = list(shape = 21), label.theme = element_text(size = 16, angle = 0, face = "italic"))) +
    #scale_fill_manual(values = c("white", "transparent")) + 
    xlim(-0.5,0.5) + 
    ylim(-0.6,0.6) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0), size = 2),
          axis.text.x = element_text(angle = 0, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = 20),
          title = element_text(size = 20)) 

#remove layers to remove the dots from plot_ordination
p$layers
p$layers <- p$layers[2]
p

#with desity lines in the plot
pdf(file = "plots/Beta_spanish_manuscript_dens.pdf", height = 6, width = 6, useDingbats = F)
plot_ordination(x_spa, 
                ordinate(x_spa, method="PCoA", "unifrac"), 
                type="sample", 
                color="Hostplant",
                shape="Hostplant", #label="ID",
                title="Catalonia", axes=c(1,2), justDF=F) + 
    facet_grid(~Hostplant) +
    geom_density2d() +
    scale_shape_manual(values = c(21,19),
                       guide=
                           guide_legend(label.theme = element_text(angle = 0, face = "italic", size = 10))) +
    scale_color_manual(values = c("orangered3", "orangered3")) +
    geom_point(aes(stroke=1, fill=Hostplant), size=4) +
    scale_fill_manual(values = c("white", "transparent")) + 
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0)),
          axis.text.x = element_text(angle = 0, size = 12),
          axis.text.y = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key = element_blank()) 
dev.off()


#change order of levels
sample_data(x_swe)$Hostplant <- factor(x=sample_data(x_swe)$Hostplant, levels = c("L.dorycnium","L.corniculatus"))
sample_data(x_swe)$Hostplant
levels(sample_data(x_swe)$Hostplant) <- c("L. dorycnium", "L. corniculatus")

#green
p_swe <-plot_ordination(x_swe, 
                    ordinate(x_swe, method="PCoA", "unifrac"), 
                    type="sample", 
                    #color="Hostplant",
                    #label="ID",
                    title="Sweden", axes=c(1,2), justDF=F#, 
                    #position=position_jitter(seed=10, width = 0.5, height = 0)
) +
    geom_point(aes(fill=Hostplant),
               size=4,
               shape=21, 
               position=position_jitter(width = 0.05, height = 0.05, seed = 10)
    ) +
    #scale_shape_manual(values = c(21,19),
    scale_fill_manual(values = c("chartreuse4", "chartreuse2"),
                      guide=
                          guide_legend(override.aes = list(shape = 21), label.theme = element_text(size = 16, angle = 0, face = "italic"))) +
    #scale_fill_manual(values = c("white", "transparent")) + 
    xlim(-0.5,0.5) + 
    ylim(-0.6,0.6) +
    theme(panel.border = element_rect(color = "black", fill = alpha("white", 0), size = 2),
          axis.text.x = element_text(angle = 0, size = 16),
          axis.text.y = element_text(size = 16),
          axis.title = element_text(size = 16),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          legend.key = element_blank(),
          legend.title = element_text(size = 20),
          title = element_text(size = 20)) 

#remove layers to remove the dots from plot_ordination
p_swe$layers
p_swe$layers <- p_swe$layers[2]
p_swe

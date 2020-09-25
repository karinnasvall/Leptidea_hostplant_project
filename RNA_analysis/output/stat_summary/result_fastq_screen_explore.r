#this script explore the conaminant in the rnaseq data, the result from FASTQ screen

setwd("Desktop/RNAseq_microbiome/")
result_fastq_screen <- read.table("mqc_fastq_screen-1_1.txt", header = T, row.names = 1)
str(result_fastq_screen)
summary(result_fastq_screen)
colSums(result_fastq_screen)
#in percentage
#round(df/rowSums(df), 2)
result_fastq_screen_percent <- round(result_fastq_screen/rowSums(result_fastq_screen),2)
#without the gutsamples (2x6 last samples)
result_fastq_screen_percent_abd <- result_fastq_screen_percent[1:130,]
summary(result_fastq_screen_percent_abd)

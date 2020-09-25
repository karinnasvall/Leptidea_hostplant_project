#table for sign GO terms
Summary_GO_MF <- read.table("Result/GOsummary/Summary_GO_MF.table", header = T)
Summary_GO_CC <- read.table("Result/GOsummary/Summary_GO_CC.table", header = T)
Summary_GO_BP <- read.table("Result/GOsummary/Summary_GO_BP.table", header = T)
Summary_GO_BP
#remove unneccessary rows and sort by genenames
Summary_GO_MF_clean <- subset(Summary_GO_MF, Summary_GO_MF$X!="NA")
Summary_GO_MF_clean_sorted <- Summary_GO_MF_clean[order(Summary_GO_MF_clean$GeneName), ]
Summary_GO_MF_clean_sorted

Summary_GO_CC_clean <- subset(Summary_GO_CC, Summary_GO_CC$X!="NA")
Summary_GO_CC_clean_sorted <- Summary_GO_CC_clean[order(Summary_GO_CC_clean$GeneName), ]
Summary_GO_CC_clean_sorted

Summary_GO_BP_clean <- subset(Summary_GO_BP, Summary_GO_BP$X!="NA")
Summary_GO_BP_clean_sorted <- Summary_GO_BP_clean[order(Summary_GO_BP_clean$GeneName), ]
Summary_GO_BP_clean_sorted


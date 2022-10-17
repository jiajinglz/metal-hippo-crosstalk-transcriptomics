###download DEseq2###
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("DESeq2")
#########################################
###import counts (integers) data and format matrix##
countsTable <- read.csv(file="C:/Users/JJ/Desktop/han_metal/counts_summary_use.csv",header=T,stringsAsFactors =FALSE)
countsMatrix <- data.matrix(countsTable[,2:ncol(countsTable)])
rownames(countsMatrix) <- countsTable[,1] #set 1st column as row names for matrix
######################
###set up experiment##
library(DESeq2)
conds <- factor(c("wt_ut","wt_ut","wt_ut","wt_cd","wt_cd","wt_cd","wt_zn","wt_zn","wt_zn",
                  "dko_ut","dko_ut","dko_ut","dko_cd","dko_cd","dko_cd","dko_zn","dko_zn","dko_zn"))
summary(conds)
colData=data.frame(condition=conds)
dds <- DESeqDataSetFromMatrix(countsMatrix,colData,design=~condition)
###variance and PCA###
rld<-rlog(dds)
colnames(rld)=colnames(countsMatrix)
plotPCA(rld, intgroup = "condition")
############################################
###################DESeq2###################
#########obtain normalized counts###########
#########and log transformed values#########
############################################
wtcounts <- read.table("C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/raw_counts/wt_counts.txt", stringsAsFactors =FALSE)
dkocounts <- read.table("C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/raw_counts/dko_counts.txt", stringsAsFactors =FALSE)
counts <- as.data.frame(matrix(nrow= nrow(wtcounts), ncol=15))
rownames(counts)= wtcounts[,1]
counts[, 1:9]= wtcounts[, 2:10]
counts[, 10:15]= dkocounts[, 5:10]
conds <- factor(c("wt", "wt", "wt", "cd", "cd", "cd", "zn", "zn", "zn",
                  "dko_cd", "dko_cd", "dko_cd", "dko_zn", "dko_zn", "dko_zn"))
library(DESeq2)
colData= data.frame(condition= conds)
dds <- DESeqDataSetFromMatrix(counts,colData,design=~condition)
sf_dds <- estimateSizeFactors(dds)
normalized_counts <- data.frame(counts(sf_dds, normalized= T), stringsAsFactors = F)
normalized_counts1 <- normalized_counts[rowMeans(normalized_counts)>2,]
ensembl_gene_id= data.frame(rownames(normalized_counts1), stringsAsFactors = FALSE)
normalized_counts2 <- cbind(ensembl_gene_id, normalized_counts1)
colnames(normalized_counts2) <- c("ensembl_gene_id", "wt_1", "wt_2", "wt_3", 
                                  "cd_1", "cd_2", "cd_3", "zn_1", "zn_2", "zn_3", 
                                  "dko_cd_1", "dko_cd_2", "dko_cd_3", 
                                  "dko_zn_1", "dko_zn_2", "dko_zn_3")
rld <- rlog(dds)
rld <- data.frame(assay(rld), stringsAsFactors = F)
ensembl_gene_id= data.frame(rownames(rld), stringsAsFactors = FALSE)
rld1 <- cbind(ensembl_gene_id, rld)
colnames(rld1) <- c("ensembl_gene_id", "wt_1", "wt_2", "wt_3", 
                    "cd_1", "cd_2", "cd_3", "zn_1", "zn_2", "zn_3", 
                    "dko_cd_1", "dko_cd_2", "dko_cd_3", 
                    "dko_zn_1", "dko_zn_2", "dko_zn_3")
###################DESeq2###################
#########DE genes for Cd and heatmap########
############################################
counts_wtvscd <- as.data.frame(matrix(nrow= nrow(counts), ncol=6))
rownames(counts_wtvscd) <- rownames(counts)
counts_wtvscd <- counts[,1:6]
library(DESeq2)
conds_wtvscd <- factor(c("wt", "wt", "wt", "cd", "cd", "cd"))
colData_wtvscd = data.frame(condition= conds_wtvscd)
dds_wtvscd <- DESeqDataSetFromMatrix(counts_wtvscd,colData_wtvscd,design=~condition)
dds_wtvscd$condition <- relevel(dds_wtvscd$condition, "wt")
dds2_wtvscd <- DESeq(dds_wtvscd)
deseq_result_wtvscd <- results(dds2_wtvscd)
####plotMA(deseq_result_wtvscd, ylim= c(-5,10))
degenes_wtvscd <- deseq_result_wtvscd[intersect(rownames(deseq_result_wtvscd[which(abs(deseq_result_wtvscd$log2FoldChange)>1),]),                                               
                                                rownames(deseq_result_wtvscd[which(deseq_result_wtvscd$padj<0.05),])),]
write.table(as.data.frame(degenes_wtvscd),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/degenes_cdvswt.txt", quote = FALSE, sep= "\t")
degenesid_wtvscd <- data.frame(paste(rownames(degenes_wtvscd)), stringsAsFactors = F)
colnames(degenesid_wtvscd) <- c("ensembl_gene_id")
cd_degene_rld <- dplyr::semi_join (rld1[,1:7], degenesid_wtvscd, by= "ensembl_gene_id")
rownames(cd_degene_rld) <- cd_degene_rld[,1]
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
cd_vs_wt_hp <- pheatmap(cd_degene_rld[,-1], 
                        scale = "row",
                        cutree_rows = 2, 
                        cluster_cols = T, 
                        show_rownames = F, 
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_cd_vs_wt <- data.frame(cutree(cd_vs_wt_hp$tree_row, k=2))
write.table(as.data.frame(genes_in_cluster_cd_vs_wt),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/genes_in_clusters_cd_vs_wt_heatmap.txt", quote = FALSE, col.names = F)
row_label_cd_vs_wt_hp <- rownames(cd_degene_rld[cd_vs_wt_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_cd_vs_wt_hp),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/row_label_cd_vs_wt_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
###################DESeq2###################
#########DE genes for Zn and heatmap########
############################################
counts_wtvszn <- as.data.frame(matrix(nrow= nrow(counts), ncol=6))
rownames(counts_wtvszn) <- rownames(counts)
counts_wtvszn[,1:3] <- counts[,1:3]
counts_wtvszn[,4:6] <- counts[,7:9]
library(DESeq2)
conds_wtvszn <- factor(c("wt", "wt", "wt", "zn", "zn", "zn"))
colData_wtvszn = data.frame(condition= conds_wtvszn)
dds_wtvszn <- DESeqDataSetFromMatrix(counts_wtvszn,colData_wtvszn,design=~condition)
dds_wtvszn$condition <- relevel(dds_wtvszn$condition, "wt")
dds2_wtvszn <- DESeq(dds_wtvszn)
deseq_result_wtvszn <- results(dds2_wtvszn)
####plotMA(deseq_result_wtvszn, ylim= c(-5,10))
degenes_wtvszn <- deseq_result_wtvszn[intersect(rownames(deseq_result_wtvszn[which(abs(deseq_result_wtvszn$log2FoldChange)>1),]),                                               
                                                rownames(deseq_result_wtvszn[which(deseq_result_wtvszn$padj<0.05),])),]
write.table(as.data.frame(degenes_wtvszn),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/degenes_znvswt.txt", quote = FALSE, sep= "\t")
degenesid_wtvszn <- data.frame(paste(rownames(degenes_wtvszn)), stringsAsFactors = F)
colnames(degenesid_wtvszn) <- c("ensembl_gene_id")
zn_degene_rld <- dplyr::semi_join (rld1[,c(1:4,8:10)], degenesid_wtvszn, by= "ensembl_gene_id")
rownames(zn_degene_rld) <- zn_degene_rld[,1]
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
zn_vs_wt_hp <- pheatmap(zn_degene_rld[,-1], 
                        scale = "row",
                        cutree_rows = 2, 
                        cluster_cols = T, 
                        show_rownames = F, 
                        clustering_method = "ward.D2", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_zn_vs_wt <- data.frame(cutree(zn_vs_wt_hp$tree_row, k=2))
write.table(as.data.frame(genes_in_cluster_zn_vs_wt),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/genes_in_clusters_zn_vs_wt_heatmap.txt", quote = FALSE, col.names = F)
row_label_zn_vs_wt_hp <- rownames(zn_degene_rld[zn_vs_wt_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_zn_vs_wt_hp),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/row_label_zn_vs_wt_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
##################################################
######find Zn and Cd overlap or pooled genes######
######extract rlog transformed value##############
##################################################
degenes_cdzn_overlap <- data.frame(intersect(degenesid_wtvscd$ensembl_gene_id, degenesid_wtvszn$ensembl_gene_id), stringsAsFactors = FALSE)
colnames(degenes_cdzn_overlap) <- c("ensembl_gene_id")
cdzn_overlap_rcounts <- dplyr::semi_join (rld1, degenes_cdzn_overlap, by= "ensembl_gene_id")
rownames(cdzn_overlap_rcounts) <- cdzn_overlap_rcounts[,1]
cdzn_overlap_rcounts <- as.matrix(cdzn_overlap_rcounts[,-1])
##################################################
#####heatmap for overlap zn and cd lists#########
##################################################
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
cdzn_overlap_hp <- pheatmap(cdzn_overlap_rcounts[,4:ncol(cdzn_overlap_rcounts)], 
                            scale = "row",
                            cutree_rows = 5, 
                            cluster_cols = F, 
                            show_rownames = F, 
                            clustering_method = "ward.D", 
                            color =colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_cdzn_overlap <- data.frame(cutree(cdzn_overlap_hp$tree_row, k=5))
genes_in_cluster_cdzn_overlap <- cbind(rownames(genes_in_cluster_cdzn_overlap), genes_in_cluster_cdzn_overlap)
write.table(as.data.frame(genes_in_cluster_cdzn_overlap),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/genes_in_clusters_cdzn_overlaplist_heatmap.txt", quote = FALSE, col.names = F, row.names= F)
row_label_cdzn_overlap_hp <- rownames(cdzn_overlap_rcounts[cdzn_overlap_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_cdzn_overlap_hp),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/row_label_cdzn_overlaplist_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
row_label_cdzn_overlap_hp_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= row_label_cdzn_overlap_hp, mart= ensembl)
write.table(as.data.frame(row_label_cdzn_overlap_hp_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/row_geneinfo_cdzn_overlaplist_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_cdzn_overlap <- unique(genes_in_cluster_cdzn_overlap[order(match(rownames(genes_in_cluster_cdzn_overlap), row_label_cdzn_overlap_hp)),])
write.table(as.data.frame(order_of_gene_cluster_cdzn_overlap),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/cdzn_overlap_heatmap_cluster_order.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_cdzn_overlap1 <- genes_in_cluster_cdzn_overlap[order(match(rownames(genes_in_cluster_cdzn_overlap), row_label_cdzn_overlap_hp)),]
extract_cdzn_overlap_cluster2_genes <- rownames(order_of_gene_cluster_cdzn_overlap1[which(order_of_gene_cluster_cdzn_overlap1[,2]==2),])
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
cdzn_overlap_cluster2_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= extract_cdzn_overlap_cluster2_genes, mart= ensembl)
write.table(as.data.frame(cdzn_overlap_cluster2_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/cdzn_overlap_cluster2_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
########################################################
######find DEgenes for Cd treated in wt vs latsdko######
########################################################
counts1 <- data.frame(cbind(rownames(counts),counts), stringsAsFactors = F)
colnames(counts1) <- c("ensembl_gene_id", "wt_1", "wt_2", "wt_3", 
                                  "cd_1", "cd_2", "cd_3", "zn_1", "zn_2", "zn_3", 
                                  "dko_cd_1", "dko_cd_2", "dko_cd_3", 
                                  "dko_zn_1", "dko_zn_2", "dko_zn_3")
DEgenes_cdzn_overlap_counts <- dplyr::semi_join (counts1, degenes_cdzn_overlap, by= "ensembl_gene_id")
cd_wtvsdko_counts <- data.frame(cbind(DEgenes_cdzn_overlap_counts[,1], DEgenes_cdzn_overlap_counts[,5:7], DEgenes_cdzn_overlap_counts[,11:13]),stringsAsFactors = F)
rownames(cd_wtvsdko_counts) <- cd_wtvsdko_counts[,1]
cd_wtvsdko_counts <- cd_wtvsdko_counts[,-1]
library(DESeq2)
conds_cd_wtvsdko<- factor(c("cdwt", "cdwt", "cdwt", "cddko", "cddko", "cddko"))
colData_cd_wtvsdko = data.frame(condition= conds_cd_wtvsdko)
dds_cd_wtvsdko <- DESeqDataSetFromMatrix(cd_wtvsdko_counts, colData_cd_wtvsdko, design=~condition)
dds_cd_wtvsdko$condition <- relevel(dds_cd_wtvsdko$condition, "cdwt")
dds2_cd_wtvsdko <- DESeq(dds_cd_wtvsdko)
deseq_result_cd_wtvsdko <- results(dds2_cd_wtvsdko)
write.table(as.data.frame(deseq_result_cd_wtvsdko),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/deseq_result_cd_wtvsdko.txt", quote = FALSE, sep= "\t")
degenes_cd_wtvsdko <- deseq_result_cd_wtvsdko[intersect(rownames(deseq_result_cd_wtvsdko[which(abs(deseq_result_cd_wtvsdko$log2FoldChange)>1),]),                                               
                                                rownames(deseq_result_cd_wtvsdko[which(deseq_result_cd_wtvsdko$padj<0.05),])),]
degenesid_cd_wtvsdko <- data.frame(paste(rownames(degenes_cd_wtvsdko)), stringsAsFactors = F)
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
cd_wtvsdko_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= degenesid_cd_wtvsdko[,1], mart= ensembl)
write.table(as.data.frame(cd_wtvsdko_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/cd_wtvdko_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
########################################################
######find DEgenes for Zn treated in wt vs latsdko######
########################################################
zn_wtvsdko_counts <- data.frame(cbind(DEgenes_cdzn_overlap_counts[,1], DEgenes_cdzn_overlap_counts[,8:10], DEgenes_cdzn_overlap_counts[,14:16]),stringsAsFactors = F)
rownames(zn_wtvsdko_counts) <- zn_wtvsdko_counts[,1]
zn_wtvsdko_counts <- zn_wtvsdko_counts[,-1]
library(DESeq2)
conds_zn_wtvsdko<- factor(c("znwt", "znwt", "znwt", "zndko", "zndko", "zndko"))
colData_zn_wtvsdko = data.frame(condition= conds_zn_wtvsdko)
dds_zn_wtvsdko <- DESeqDataSetFromMatrix(zn_wtvsdko_counts, colData_zn_wtvsdko, design=~condition)
dds_zn_wtvsdko$condition <- relevel(dds_zn_wtvsdko$condition, "znwt")
dds2_zn_wtvsdko <- DESeq(dds_zn_wtvsdko)
deseq_result_zn_wtvsdko <- results(dds2_zn_wtvsdko)
write.table(as.data.frame(deseq_result_zn_wtvsdko),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/deseq_result_zn_wtvsdko.txt", quote = FALSE, sep= "\t")
degenes_zn_wtvsdko <- deseq_result_zn_wtvsdko[intersect(rownames(deseq_result_zn_wtvsdko[which(abs(deseq_result_zn_wtvsdko$log2FoldChange)>1),]),                                               
                                                        rownames(deseq_result_zn_wtvsdko[which(deseq_result_zn_wtvsdko$padj<0.05),])),]
degenesid_zn_wtvsdko <- data.frame(paste(rownames(degenes_zn_wtvsdko)), stringsAsFactors = F)
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
zn_wtvsdko_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= degenesid_zn_wtvsdko[,1], mart= ensembl)
write.table(as.data.frame(zn_wtvsdko_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/zn_wtvdko_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
#######################################################################
######find DEgenes regulated by lats1/2 in both Cd and Zn treated######
################################Heatmap################################
double_degenes <- data.frame(intersect(degenesid_zn_wtvsdko[,1], degenesid_cd_wtvsdko[,1]), stringsAsFactors = FALSE)
colnames(double_degenes) <- c("ensembl_gene_id")
double_degenes_rcounts <- dplyr::semi_join (rld1, double_degenes, by= "ensembl_gene_id")
rownames(double_degenes_rcounts) <- double_degenes_rcounts[,1]
double_degenes_rcounts <- as.matrix(double_degenes_rcounts[,-1])
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
double_degenes_hp <- pheatmap(double_degenes_rcounts[,4:ncol(double_degenes_rcounts)], 
                            scale = "row",
                            cutree_rows = 2, 
                            cluster_cols = F, 
                            show_rownames = F, 
                            clustering_method = "ward.D", 
                            color =colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_double_degenes <- data.frame(cutree(double_degenes_hp$tree_row, k=2))
genes_in_cluster_double_degenes <- cbind(rownames(genes_in_cluster_double_degenes), genes_in_cluster_double_degenes)
write.table(as.data.frame(genes_in_cluster_double_degenes),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/genes_in_cluster_double_degenes.txt", quote = FALSE, col.names = F, row.names= F)
row_label_double_degenes_hp <- rownames(double_degenes_rcounts[double_degenes_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_double_degenes_hp),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/row_label_double_degenes_hp.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
row_label_double_degenes_hp_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= row_label_double_degenes_hp, mart= ensembl)
write.table(as.data.frame(row_label_double_degenes_hp_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/row_label_double_degenes_hp_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_double_degenes <- unique(genes_in_cluster_double_degenes[order(match(rownames(genes_in_cluster_double_degenes), row_label_double_degenes_hp)),])
write.table(as.data.frame(order_of_gene_cluster_double_degenes),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/order_of_gene_cluster_double_degenes.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
extract_double_degenes_cluster1_genes <- rownames(order_of_gene_cluster_double_degenes[which(order_of_gene_cluster_double_degenes[,2]==1),])
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
double_degenes_cluster1_downregulated_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= extract_double_degenes_cluster1_genes, mart= ensembl)
write.table(as.data.frame(double_degenes_cluster1_downregulated_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/double_degenes_cluster1_downregulated_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
extract_double_degenes_cluster2_genes <- rownames(order_of_gene_cluster_double_degenes[which(order_of_gene_cluster_double_degenes[,2]==2),])
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
double_degenes_cluster2_upregulated_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= extract_double_degenes_cluster2_genes, mart= ensembl)
write.table(as.data.frame(double_degenes_cluster2_upregulated_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/doubleDE/double_degenes_cluster2_upregulated_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")










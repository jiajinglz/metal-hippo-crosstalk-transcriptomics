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
degenes_cdzn_combine <- data.frame(unique(c(degenesid_wtvscd$ensembl_gene_id,degenesid_wtvszn$ensembl_gene_id)), stringsAsFactors = FALSE)
colnames(degenes_cdzn_combine) <- c("ensembl_gene_id")
cdzn_combine_rcounts <- dplyr::semi_join (rld1, degenes_cdzn_combine, by= "ensembl_gene_id")
rownames(cdzn_combine_rcounts) <- cdzn_combine_rcounts[,1]
cdzn_combine_rcounts <- as.matrix(cdzn_combine_rcounts[,-1])
degenes_cdzn_overlap <- data.frame(intersect(degenesid_wtvscd$ensembl_gene_id, degenesid_wtvszn$ensembl_gene_id), stringsAsFactors = FALSE)
colnames(degenes_cdzn_overlap) <- c("ensembl_gene_id")
cdzn_overlap_rcounts <- dplyr::semi_join (rld1, degenes_cdzn_overlap, by= "ensembl_gene_id")
rownames(cdzn_overlap_rcounts) <- cdzn_overlap_rcounts[,1]
cdzn_overlap_rcounts <- as.matrix(cdzn_overlap_rcounts[,-1])
##################################################
#####heatmap for combined zn and cd lists#########
##################################################
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
cdzn_combined_hp <- pheatmap(cdzn_combine_rcounts, 
                        scale = "row",
                        cutree_rows = 4, 
                        cluster_cols = F, 
                        show_rownames = F, 
                        clustering_method = "ward.D", 
                        col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_cdzn_combine <- data.frame(cutree(cdzn_combined_hp$tree_row, k=4))
genes_in_cluster_cdzn_combine <- cbind(rownames(genes_in_cluster_cdzn_combine), genes_in_cluster_cdzn_combine)
write.table(as.data.frame(genes_in_cluster_cdzn_combine),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/genes_in_clusters_cdzn_combinedlist_heatmap.txt", quote = FALSE, col.names = F)
row_label_cdzn_combined_hp <- rownames(cdzn_combine_rcounts[cdzn_combined_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_cdzn_combined_hp),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/row_label_cdzn_combinedlist_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_cdzn_combine <- unique(genes_in_cluster_cdzn_combine[order(match(rownames(genes_in_cluster_cdzn_combine), row_label_cdzn_combined_hp)),])
write.table(as.data.frame(order_of_gene_cluster_cdzn_combine),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/cdzn_combine_heatmap_cluster_order.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_cdzn_combine1 <- genes_in_cluster_cdzn_combine[order(match(rownames(genes_in_cluster_cdzn_combine), row_label_cdzn_combined_hp)),]
extract_cdzn_combine_cluster3_genes <- rownames(order_of_gene_cluster_cdzn_combine1[which(order_of_gene_cluster_cdzn_combine1[,2]==3),])
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
cdzn_combine_cluster3_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= extract_cdzn_combine_cluster3_genes, mart= ensembl)
write.table(as.data.frame(cdzn_combine_cluster3_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/cdzn_combine_cluster3_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
##################################################
#####heatmap for overlap zn and cd lists#########
##################################################
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
cdzn_overlap_hp <- pheatmap(cdzn_overlap_rcounts, 
                             scale = "row",
                             cutree_rows = 6, 
                             cluster_cols = F, 
                             show_rownames = F, 
                             clustering_method = "ward.D", 
                             color =colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_cdzn_overlap <- data.frame(cutree(cdzn_overlap_hp$tree_row, k=6))
genes_in_cluster_cdzn_overlap <- cbind(rownames(genes_in_cluster_cdzn_overlap), genes_in_cluster_cdzn_overlap)
write.table(as.data.frame(genes_in_cluster_cdzn_overlap),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/genes_in_clusters_cdzn_overlaplist_heatmap.txt", quote = FALSE, col.names = F)
row_label_cdzn_overlap_hp <- rownames(cdzn_overlap_rcounts[cdzn_overlap_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_cdzn_overlap_hp),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/row_label_cdzn_overlaplist_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_cdzn_overlap <- unique(genes_in_cluster_cdzn_overlap[order(match(rownames(genes_in_cluster_cdzn_overlap), row_label_cdzn_overlap_hp)),])
write.table(as.data.frame(order_of_gene_cluster_cdzn_overlap),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/cdzn_overlap_heatmap_cluster_order.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
order_of_gene_cluster_cdzn_overlap1 <- genes_in_cluster_cdzn_overlap[order(match(rownames(genes_in_cluster_cdzn_overlap), row_label_cdzn_overlap_hp)),]
extract_cdzn_overlap_cluster2_genes <- rownames(order_of_gene_cluster_cdzn_overlap1[which(order_of_gene_cluster_cdzn_overlap1[,2]==2),])
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
cdzn_overlap_cluster2_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= extract_cdzn_overlap_cluster2_genes, mart= ensembl)
write.table(as.data.frame(cdzn_overlap_cluster2_geneinfo),file="C:/Users/JJ/Desktop/Lab/Bioinformatics/han_metal/Differential gene expression/cdzn_overlap_cluster2_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")












############################################
###################DESeq2###################
#########obtain normalized counts###########
#########and log transformed values#########
############################################
raw_counts_table <- read.delim("C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/raw_counts/raw_counts_table.txt", header=FALSE)
counts <- as.data.frame(raw_counts_table)
rownames(counts)= counts[,1]
counts <- counts[,-1]
conds <- factor(c("wt", "wt", "wt", "cd", "cd", "cd", "zn", "zn", "zn",
                  "dko_ut", "dko_ut", "dko_ut", "dko_cd", "dko_cd", "dko_cd", 
                  "dko_zn", "dko_zn", "dko_zn", "tko_ut", "tko_ut", "tko_ut", 
                  "tko_cd", "tko_cd", "tko_cd", "tko_zn", "tko_zn", "tko_zn"))
library(DESeq2)
colData= data.frame(condition= conds)
dds <- DESeqDataSetFromMatrix(counts,colData,design=~condition)
sf_dds <- estimateSizeFactors(dds)
normalized_counts <- data.frame(counts(sf_dds, normalized= T), stringsAsFactors = F)
write.table(as.data.frame(normalized_counts),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/normalized_counts_wo_filter.txt", quote = FALSE, row.names= T, col.names= conds, sep= "\t")
rld <- rlog(dds)
plotPCA(rld, intgroup = "condition")
rld <- data.frame(assay(rld), stringsAsFactors = F)
ensembl_gene_id= data.frame(rownames(rld), stringsAsFactors = FALSE)
rld1 <- cbind(ensembl_gene_id, rld)
colnames(rld1) <- c("ensembl_gene_id", "wt_1", "wt_2", "wt_3", 
                    "cd_1", "cd_2", "cd_3", "zn_1", "zn_2", "zn_3", "dko_ut_1", "dko_ut_2",
                    "dko_ut_3", "dko_cd_1", "dko_cd_2", "dko_cd_3", "dko_zn_1", "dko_zn_2", 
                    "dko_zn_3", "tko_ut_1", "tko_ut_2", "tko_ut_3", "tko_cd_1", "tko_cd_2", 
                    "tko_cd_3", "tko_zn_1", "tko_zn_2", "tko_zn_3")
write.table(as.data.frame(rld1[,2:ncol(rld1)]),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/normalized_counts_log_transformed.txt", quote = FALSE, row.names= T, col.names= conds, sep= "\t")
###################DESeq2###################
#########DE genes for WT: UT versus Cd######
############################################
counts_wtvscd <- as.data.frame(matrix(nrow= nrow(counts), ncol=6))
rownames(counts_wtvscd) <- rownames(counts)
wtvscdCol <- c(1,2,3,4,5,6)
counts_wtvscd <- counts[,wtvscdCol]
library(DESeq2)
conds_wtvscd <- factor(c("wt", "wt", "wt", "cd", "cd", "cd"))
colData_wtvscd = data.frame(condition= conds_wtvscd)
dds_wtvscd<- DESeqDataSetFromMatrix(counts_wtvscd,colData_wtvscd,design=~condition)
dds_wtvscd$condition <- relevel(dds_wtvscd$condition, "wt")
dds2_wtvscd <- DESeq(dds_wtvscd)
deseq_result_wtvscd <- results(dds2_wtvscd)
plotMA(deseq_result_wtvscd, ylim= c(-10,10))
degenes_wtvscd <- deseq_result_wtvscd[intersect(rownames(deseq_result_wtvscd[which(abs(deseq_result_wtvscd$log2FoldChange)>1),]),                                               
                                                rownames(deseq_result_wtvscd[which(deseq_result_wtvscd$padj<0.05),])),]
write.table(as.data.frame(degenes_wtvscd),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/degenes_cdvswt.txt", quote = FALSE, sep= "\t")
degenesid_wtvscd <- data.frame(paste(rownames(degenes_wtvscd)), stringsAsFactors = F)
colnames(degenesid_wtvscd) <- c("ensembl_gene_id")
wtvscdRLDcol <- c(1,2,3,4,5,6,7)
cd_degene_rld <- dplyr::semi_join (rld1[,wtvscdRLDcol], degenesid_wtvscd, by= "ensembl_gene_id")
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
write.table(as.data.frame(genes_in_cluster_cd_vs_wt),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/genes_in_clusters_cd_vs_wt_heatmap.txt", quote = FALSE, col.names = F)
row_label_cd_vs_wt_hp <- rownames(cd_degene_rld[cd_vs_wt_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_cd_vs_wt_hp),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/row_label_cd_vs_wt_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
###################DESeq2###################
#########DE genes for WT: UT versus Zn######
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
plotMA(deseq_result_wtvszn, ylim= c(-10,10))
degenes_wtvszn <- deseq_result_wtvszn[intersect(rownames(deseq_result_wtvszn[which(abs(deseq_result_wtvszn$log2FoldChange)>1),]),                                               
                                                rownames(deseq_result_wtvszn[which(deseq_result_wtvszn$padj<0.05),])),]
write.table(as.data.frame(degenes_wtvszn),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/degenes_znvswt.txt", quote = FALSE, sep= "\t")
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
write.table(as.data.frame(genes_in_cluster_zn_vs_wt),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/genes_in_clusters_zn_vs_wt_heatmap.txt", quote = FALSE, col.names = F)
row_label_zn_vs_wt_hp <- rownames(zn_degene_rld[zn_vs_wt_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_zn_vs_wt_hp),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/row_label_zn_vs_wt_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
##################################################
######find Zn and Cd overlap or pooled genes######
######extract rlog transformed values#############
##################################################
degenes_cdzn_overlap <- data.frame(intersect(degenesid_wtvscd$ensembl_gene_id, degenesid_wtvszn$ensembl_gene_id), stringsAsFactors = FALSE)
colnames(degenes_cdzn_overlap) <- c("ensembl_gene_id")
cdzn_overlap_rld <- dplyr::semi_join (rld1, degenes_cdzn_overlap, by= "ensembl_gene_id")
rownames(cdzn_overlap_rld) <- cdzn_overlap_rld[,1]
cdzn_overlap_rld <- as.matrix(cdzn_overlap_rld[,-1])
##########################################################
#####heatmap for metal (Cd, Zn) responsisve genes#########
##########################################################
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
cdzn_overlap_hp <- pheatmap(cdzn_overlap_rld, 
                            scale = "row",
                            cutree_rows = 3, 
                            cluster_cols = F, 
                            show_rownames = F, 
                            clustering_method = "ward.D", 
                            color =colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
colwtvstko <- c(1:9, 19:27)
cdzn_overlap_hp_wtvsTko <- pheatmap(cdzn_overlap_rld[,colwtvstko], 
                            scale = "row",
                            cutree_rows = 3, 
                            cluster_cols = F, 
                            show_rownames = F, 
                            clustering_method = "ward.D", 
                            color =colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
genes_in_cluster_cdzn_overlap <- data.frame(cutree(cdzn_overlap_hp$tree_row, k=3))
genes_in_cluster_cdzn_overlap <- cbind(rownames(genes_in_cluster_cdzn_overlap), genes_in_cluster_cdzn_overlap)
write.table(as.data.frame(genes_in_cluster_cdzn_overlap),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/genes_in_clusters_cdzn_overlaplist_heatmap.txt", quote = FALSE, col.names = F, row.names= F)
row_label_cdzn_overlap_hp <- rownames(cdzn_overlap_rld[cdzn_overlap_hp$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
write.table(as.data.frame(row_label_cdzn_overlap_hp),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/row_label_cdzn_overlaplist_heatmap.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
#####I skip annotation for genes in cluster part######
######scripts available in double DE floder###########
########################################################
####find metal responsive genes changed upon latsdko####
#################when treated with Cd###################
########################################################
counts1 <- data.frame(cbind(rownames(counts),counts), stringsAsFactors = F)
colnames(counts1) <- c("ensembl_gene_id", "wt_1", "wt_2", "wt_3", 
                       "cd_1", "cd_2", "cd_3", "zn_1", "zn_2", "zn_3", "dko_ut_1", "dko_ut_2",
                       "dko_ut_3", "dko_cd_1", "dko_cd_2", "dko_cd_3", "dko_zn_1", "dko_zn_2", 
                       "dko_zn_3", "tko_ut_1", "tko_ut_2", "tko_ut_3", "tko_cd_1", "tko_cd_2", 
                       "tko_cd_3", "tko_zn_1", "tko_zn_2", "tko_zn_3")
DEgenes_cdzn_overlap_counts <- dplyr::semi_join (counts1, degenes_cdzn_overlap, by= "ensembl_gene_id")
cd_wtvsdko_counts <- data.frame(cbind(DEgenes_cdzn_overlap_counts[,1], DEgenes_cdzn_overlap_counts[,5:7], DEgenes_cdzn_overlap_counts[,14:16]),stringsAsFactors = F)
rownames(cd_wtvsdko_counts) <- cd_wtvsdko_counts[,1]
cd_wtvsdko_counts <- cd_wtvsdko_counts[,-1]
library(DESeq2)
conds_cd_wtvsdko<- factor(c("cdwt", "cdwt", "cdwt", "cddko", "cddko", "cddko"))
colData_cd_wtvsdko = data.frame(condition= conds_cd_wtvsdko)
dds_cd_wtvsdko <- DESeqDataSetFromMatrix(cd_wtvsdko_counts, colData_cd_wtvsdko, design=~condition)
sf_dds_cd_wtvsdko <- estimateSizeFactors(dds_cd_wtvsdko)
normalized_counts_cd_wtvsdko <- data.frame(counts(sf_dds_cd_wtvsdko, normalized= T), stringsAsFactors = F)
###looks like this normalized counts with just metal responsive genes are different from all genes#####
####can test DE genes bewteen cd wt and dko for all genes to see if overlap with result from###########
#########using just the metal responsive genes#########################################################
dds_cd_wtvsdko$condition <- relevel(dds_cd_wtvsdko$condition, "cdwt")
dds2_cd_wtvsdko <- DESeq(dds_cd_wtvsdko)
deseq_result_cd_wtvsdko <- results(dds2_cd_wtvsdko)
write.table(as.data.frame(deseq_result_cd_wtvsdko),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/deseq_result_cd_wtvsdko.txt", quote = FALSE, sep= "\t")
degenes_cd_wtvsdko <- deseq_result_cd_wtvsdko[intersect(rownames(deseq_result_cd_wtvsdko[which(abs(deseq_result_cd_wtvsdko$log2FoldChange)>1),]),                                               
                                                        rownames(deseq_result_cd_wtvsdko[which(deseq_result_cd_wtvsdko$padj<0.05),])),]
degenesid_cd_wtvsdko <- data.frame(paste(rownames(degenes_cd_wtvsdko)), stringsAsFactors = F)
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
cd_wtvsdko_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= degenesid_cd_wtvsdko[,1], mart= ensembl)
write.table(as.data.frame(cd_wtvsdko_geneinfo),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/cd_wtvdko_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
#############################################################
######find metal responsive genes regulated by lats1/2#######
####################when treated with Zn#####################
zn_wtvsdko_counts <- data.frame(cbind(DEgenes_cdzn_overlap_counts[,1], DEgenes_cdzn_overlap_counts[,8:10], DEgenes_cdzn_overlap_counts[,17:19]),stringsAsFactors = F)
rownames(zn_wtvsdko_counts) <- zn_wtvsdko_counts[,1]
zn_wtvsdko_counts <- zn_wtvsdko_counts[,-1]
library(DESeq2)
conds_zn_wtvsdko<- factor(c("znwt", "znwt", "znwt", "zndko", "zndko", "zndko"))
colData_zn_wtvsdko = data.frame(condition= conds_zn_wtvsdko)
dds_zn_wtvsdko <- DESeqDataSetFromMatrix(zn_wtvsdko_counts, colData_zn_wtvsdko, design=~condition)
dds_zn_wtvsdko$condition <- relevel(dds_zn_wtvsdko$condition, "znwt")
dds2_zn_wtvsdko <- DESeq(dds_zn_wtvsdko)
deseq_result_zn_wtvsdko <- results(dds2_zn_wtvsdko)
write.table(as.data.frame(deseq_result_zn_wtvsdko),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/deseq_result_zn_wtvsdko.txt", quote = FALSE, sep= "\t")
degenes_zn_wtvsdko <- deseq_result_zn_wtvsdko[intersect(rownames(deseq_result_zn_wtvsdko[which(abs(deseq_result_zn_wtvsdko$log2FoldChange)>1),]),                                               
                                                        rownames(deseq_result_zn_wtvsdko[which(deseq_result_zn_wtvsdko$padj<0.05),])),]
degenesid_zn_wtvsdko <- data.frame(paste(rownames(degenes_zn_wtvsdko)), stringsAsFactors = F)
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
zn_wtvsdko_geneinfo <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= degenesid_zn_wtvsdko[,1], mart= ensembl)
write.table(as.data.frame(zn_wtvsdko_geneinfo),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/zn_wtvdko_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
#######################################################################
######find lats1/2 reguated metal(Cd&Zn) responsive genes #############
################################Heatmap################################
double_degenes <- data.frame(intersect(degenesid_zn_wtvsdko[,1], degenesid_cd_wtvsdko[,1]), stringsAsFactors = FALSE)
colnames(double_degenes) <- c("ensembl_gene_id")
double_degenes_rld <- dplyr::semi_join (rld1, double_degenes, by= "ensembl_gene_id")
rownames(double_degenes_rld) <- double_degenes_rld[,1]
double_degenes_rld <- as.matrix(double_degenes_rld[,-1])
library(pheatmap)
library("RColorBrewer")
ramp<-1:3/3 
cols<-c(rgb(ramp,0,0),rgb(0,ramp,0),rgb(0,0,ramp),rgb(ramp,0,ramp)) 
double_degenes_hp_all <- pheatmap(double_degenes_rld, 
                              scale = "row",
                              cutree_rows = 5, 
                              cluster_cols = T, 
                              show_rownames = F, 
                              clustering_method = "ward.D", 
                              color =colorRampPalette(rev(brewer.pal(9,"RdBu")))(255))
##############################################################
###############Annotate genes in certain cluster##############
##############################################################
genes_in_cluster_double_degenes <- data.frame(cutree(double_degenes_hp_all$tree_row, k=5))
genes_in_cluster_double_degenes <- cbind(rownames(genes_in_cluster_double_degenes), genes_in_cluster_double_degenes)
row_label_double_degenes_hp_all <- rownames(double_degenes_rld[double_degenes_hp_all$tree_row[["order"]],]) #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
rowlabel_with_cluster_double_degenes_hp_all <- genes_in_cluster_double_degenes[match(row_label_double_degenes_hp_all, genes_in_cluster_double_degenes[,1]),]
extract_double_degenes_hp_all_cluster1_genes <- rowlabel_with_cluster_double_degenes_hp_all[which(rowlabel_with_cluster_double_degenes_hp_all[,2]==1),]
extract_double_degenes_hp_all_cluster2_genes <- rowlabel_with_cluster_double_degenes_hp_all[which(rowlabel_with_cluster_double_degenes_hp_all[,2]==2),]
extract_double_degenes_hp_all_cluster3_genes <- rowlabel_with_cluster_double_degenes_hp_all[which(rowlabel_with_cluster_double_degenes_hp_all[,2]==3),]
extract_double_degenes_hp_all_cluster4_genes <- rowlabel_with_cluster_double_degenes_hp_all[which(rowlabel_with_cluster_double_degenes_hp_all[,2]==4),]
extract_double_degenes_hp_all_cluster5_genes <- rowlabel_with_cluster_double_degenes_hp_all[which(rowlabel_with_cluster_double_degenes_hp_all[,2]==5),]
library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
double_degenes_hp_all_cluster2_genes_info <- getBM(filters = 'ensembl_gene_id', attributes=c('ensembl_gene_id', 'external_gene_name', 'gene_biotype'), values= extract_double_degenes_hp_all_cluster2_genes[,1], mart= ensembl)
write.table(as.data.frame(double_degenes_hp_all_cluster2_genes_info),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/double_degenes_hp_all_cluster2_genes_info.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
################################################################
####check cluster 2 metal included cluster DEseq stats##########
##################in metal responsive filter####################
double_degenes_cluster2_geneid <- as.character(extract_double_degenes_hp_all_cluster2_genes[,1])
deseq_result_wtvscd1 <- as.data.frame(cbind(rownames(deseq_result_wtvscd), deseq_result_wtvscd))
cd_de_stats_for_cluster2_genes <- deseq_result_wtvscd1[match(double_degenes_cluster2_geneid, deseq_result_wtvscd1[,1]),]
cd_de_stats_for_cluster2_genes1 <- cd_de_stats_for_cluster2_genes[match(double_degenes_hp_all_cluster2_genes_info[,1], cd_de_stats_for_cluster2_genes[,1]),]
cd_de_stats_for_cluster2_genes_with_geneinfo <- cbind(double_degenes_cdznwtdkotko_cluster2_genes_info,cd_de_stats_for_cluster2_genes1)
write.table(as.data.frame(cd_de_stats_for_cluster2_genes_with_geneinfo),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/cd_de_stats_for_cluster2_genes_with_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")
deseq_result_wtvszn1 <- as.data.frame(cbind(rownames(deseq_result_wtvszn), deseq_result_wtvszn))
zn_de_stats_for_cluster2_genes <- deseq_result_wtvszn1[match(double_degenes_cluster2_geneid, deseq_result_wtvszn1[,1]),]
zn_de_stats_for_cluster2_genes1 <- zn_de_stats_for_cluster2_genes[match(double_degenes_hp_all_cluster2_genes_info[,1], zn_de_stats_for_cluster2_genes[,1]),]
zn_de_stats_for_cluster2_genes_with_geneinfo <- cbind(double_degenes_cdznwtdkotko_cluster2_genes_info,zn_de_stats_for_cluster2_genes1)
write.table(as.data.frame(zn_de_stats_for_cluster2_genes_with_geneinfo),file="C:/Users/latro/Desktop/Lab/Bioinformatics/han_metal/DE gene_with TKO samples/zn_de_stats_for_cluster2_genes_with_geneinfo.txt", quote = FALSE, row.names= F, col.names= F, sep= "\t")

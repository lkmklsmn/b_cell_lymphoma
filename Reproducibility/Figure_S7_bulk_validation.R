setwd("/Users/yanfangfang/Downloads/MW/")
# Load R libs ####
library(Seurat)
library(dplyr)
library(nlme)
library(ggplot2)
library(DESeq2)
load("validate.RData")
res_tmp <- res[res$adj.p<0.05,]
res_tmp <- res_tmp[order(res_tmp$adj.p,decreasing = F),]
mixed_degs <- rownames(res_tmp)
Idents(seu) <- seu$ibrutinib_sensitivity
stan_df <- FindAllMarkers(seu,only.pos = T)
stan_degs <- unique(stan_df$gene)
seu <- FindVariableFeatures(seu)
high <- VariableFeatures(seu)

# Figure S4A ####
bulk <- read.delim("data/Zhangliang_STM_paper_processed_data.txt")
sens <- read.delim("data/Zhangliang_STM_paper_sensitivity.txt")
sens$sensitivity[sens$sample %in% c("MCL18","MCL20","MCL3","MCL12","MCL4","MCL2")] <- "PR"
sens <- sens[sens$sensitivity!="PR",]
bulk <- bulk[,colnames(bulk) %in% sens$sample]
# Project into PCA defined by single cell data
project_deg <- function(degs){
        expr <- seu@assays$RNA@counts[degs,]
        asplit <- split(rownames(meta), meta$sample)
        pseudobulks <- do.call(cbind, lapply(asplit, function(x) rowMeans(expr[,x])))
        # match genes
        ok <- intersect(rownames(pseudobulks),rownames(bulk))
        pseudobulks_ok <- pseudobulks[ok,]
        bulk_ok <- bulk[ok,]
        # voom 
        pseudobulks_ok <- voom(pseudobulks_ok)$E
        bulk_ok <- voom(bulk_ok)$E
        merge <- cbind(pseudobulks_ok,bulk_ok)
        library('preprocessCore')
        norm_merge <- normalize.quantiles(merge)
        colnames(norm_merge) <- colnames(merge)
        rownames(norm_merge) <- rownames(merge)
        
        # pca ####
        tmp <- norm_merge[, colnames(pseudobulks)]
        pca <- prcomp(t(tmp))
        
        # project bulk into pca sapce ####
        preds <- predict(pca, t(norm_merge[, colnames(bulk)]))
        df <- as.data.frame(rbind(pca$x,preds))
        df$group <- "pseudobulk"
        df[colnames(bulk),'group'] <- "bulk"
        df$sensitivity <- "S"
        r_sample <- c(unique(seu$sample[seu$ibrutinib_sensitivity=="R"]),sens$sample[sens$sensitivity=="R"])
        df[r_sample,'sensitivity'] <- "R"
        ggplot(df,aes(PC1,PC2,color=sensitivity,shape=group))+geom_point()+theme_bw()+
                facet_wrap(~group,scales='free')
}
p1 <- project_deg(mixed_degs)+ggtitle('PCA based on mixed model DEGs')+theme(legend.position = "none")
p2 <- project_deg(stan_degs)+ggtitle('PCA based on standard DEGs')
p3 <- project_deg(high)+ggtitle('PCA based on highly variable genes')+theme(legend.position = "none")
p1/p2/p3

# Volcano plot ####
condition <- factor(sens$sensitivity,levels=c("S","R"))
colData <- data.frame(row.names = colnames(bulk), condition)
countData <- bulk[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ condition)
dds <- DESeq(dds)
res_bulk <- results(dds)
res_bulk <- data.frame(res_bulk)

res_bulk$bulk_sig <- NA
res_bulk$bulk_sig[res_bulk$padj<0.05 & res_bulk$log2FoldChange>0] <- "bulk, up"
res_bulk$bulk_sig[res_bulk$padj<0.05 & res_bulk$log2FoldChange<0] <- "bulk, down"

res_sc <- res
res_sc$sc_mixed_sig <- NA
res_sc[rownames(res_sc)[res_sc$adj.p<0.05 & res_sc$coef_S>0],'sc_mixed_sig'] <- "sc_mixed, down"
res_sc[rownames(res_sc)[res_sc$adj.p<0.05 & res_sc$coef_S<0],'sc_mixed_sig'] <- "sc_mixed, up"

stan_up <- rownames(stan_df)[stan_df$p_val_adj<0.05 & stan_df$cluster=="R"]
stan_down <- rownames(stan_df)[stan_df$p_val_adj<0.05 & stan_df$cluster=="S"]

res_bulk$sc_mixed_sig <- NA
res_bulk[rownames(res_bulk) %in% rownames(res_sc)[res_sc$sc_mixed_sig=='sc_mixed, down'],'sc_mixed_sig'] <- "sc_mixed, down"
res_bulk[rownames(res_bulk) %in% rownames(res_sc)[res_sc$sc_mixed_sig=='sc_mixed, up'],'sc_mixed_sig'] <- "sc_mixed, up"

res_bulk$sc_standard_sig <- NA
res_bulk[rownames(res_bulk) %in% stan_up,'sc_standard_sig'] <- 'sc_standard,up'
res_bulk[rownames(res_bulk) %in% stan_down,'sc_standard_sig'] <- 'sc_standard,down'

table(res_bulk$bulk_sig)
table(res_bulk$sc_sig)
table(res_bulk$bulk_sig,res_bulk$sc_mixed_sig)
fisher.test(table(res_bulk$bulk_sig,res_bulk$sc_mixed_sig))
table(res_bulk$bulk_sig,res_bulk$sc_standard_sig)
fisher.test(table(res_bulk$bulk_sig,res_bulk$sc_standard_sig))




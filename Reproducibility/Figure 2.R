setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
library(pheatmap)
load("/Users/yanfangfang/BlueYellowColormaps_V1.RData")
integrated <- readRDS(file="data/integrated.rds")
integrated$celltype[integrated$celltype=="Tumor B" & integrated$sample %in% c("Normal_1","Normal_2")] <- "Normal B"
b_seu <- readRDS("data/integrated_b_cells.rds")
good <- names(table(b_seu$sample))[(table(b_seu$sample)>20)]
b_seu <- b_seu[,b_seu$sample %in% good]

# Figure 2A ####
integrated$tmp <- "non-B"
integrated$tmp[integrated$celltype=="Tumor B"] <- "Tumor-B"
p1 <- DimPlot(integrated,group.by = 'tmp',cols=c("lightgrey","steelblue"))+
        labs(title="B cells")+NoLegend()
b_seu$clinical_outcome <- factor(b_seu$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(b_seu$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
p2 <- DimPlot(b_seu,group.by = 'clinical_outcome',
        cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Clinical outcome")+NoLegend()
p3 <- DimPlot(b_seu,group.by = 'patient',label=T,repel=T)+
        labs(title="Patient")+NoLegend()
p4 <- DimPlot(b_seu,group.by = 'sample',label=T,repel=T)+
        labs(title="Sample")+NoLegend()
p1|p2|p3|p4

# Figure 2C ####
# infercnv_obj <- readRDS('inferCNV/run.final.infercnv_obj_allcells')
# seu <- CreateSeuratObject(infercnv_obj@expr.data)
# rm(infercnv_obj)
# meta <- read.delim("meta/meta_b_cells.txt")
# seu@meta.data <- meta[colnames(seu),]
# seu <- NormalizeData(seu)
# seu <- ScaleData(seu, features = rownames(seu))
# seu <- RunPCA(seu,features=rownames(seu))
# seu <- FindNeighbors(seu, dims = 1:10)
# seu <- FindClusters(seu, resolution = 0.4)
# seu <- RunUMAP(seu, dims = 1:10, check_duplicates = FALSE)
# seu$clinical_outcome <- factor(seu$clinical_outcomes,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
# require(graphics)
# dist.mtx <- as.matrix(dist(seu@reductions$umap@cell.embeddings,method = "euclidean",upper=TRUE))
# df <- seu@meta.data
# lst <- levels(seu$clinical_outcome)
# res <- do.call(rbind,lapply(lst,function(x) {
#         c1 <- rownames(df)[df$clinical_outcome=="Normal"]
#         c2 <- sample(rownames(df)[df$clinical_outcome==x],length(c1))
#         df <- dist.mtx[c1,c2]
#         dist <- apply(df,2,mean)
#         tmp <- data.frame('dist'=dist,x)
#         colnames(tmp) <- c("dist",'clinical_outcome')
#         tmp
# }))
# res <- data.frame(res)
# res$sample <- df[match(rownames(res),rownames(df)),'sample']
# res$source <- df[match(res$sample,df$sample),'source']
# library(ggpubr)
# res$clinical_outcome <- factor(res$clinical_outcome,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
# levels(res$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
# saveRDS(res,"data/genome_instability_res.rds")
res <- readRDS("data/genome_instability_res.rds")
ggplot(res[res$source=="PB",], aes(x=clinical_outcome,y=dist,fill=clinical_outcome)) + 
        geom_boxplot(outlier.shape=NA)+
        xlab("clincial outcome")+ylab("Euclidean distance")+
        ggtitle("Genome instability score")+
        theme_bw()+scale_fill_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        stat_compare_means(method = 'anova')

# Figure 2E, chr12,chr17 mutated genes ####
infer <- readRDS("inferCNV/run.final.infercnv_obj_allcells")
expr <- infer@expr.data
# Check chr 12 ####
geneinfo <- (infer@gene_order)
ok <- which(geneinfo$chr == "chr12")
aframe <- data.frame(expr[ok,])
means <- do.call(cbind, lapply(infer@observation_grouped_cell_indices, function(x) rowMeans(expr[ok, x])))
means <- data.frame(means)
means$position <- geneinfo$start[ok]
tmp <- reshape2::melt(means, id.vars = "position")
ggplot(tmp, aes(position, value, group = variable, color = variable)) +
        geom_point() +
        geom_hline(yintercept = 1) +
        geom_vline(xintercept = 7500000, color = "red") +
        theme_classic()
hits <- rownames(geneinfo)[geneinfo$start < 7500000 & geneinfo$chr == "chr12"]

# Check chr 17 ####
geneinfo <- (infer@gene_order)
ok <- which(geneinfo$chr == "chr17")
aframe <- data.frame(expr[ok,])
means <- do.call(cbind, lapply(infer@observation_grouped_cell_indices, function(x) rowMeans(expr[ok, x])))
means <- data.frame(means)
means$position <- geneinfo$start[ok]
tmp <- reshape2::melt(means, id.vars = "position")
ggplot(tmp, aes(position, value, group = variable, color = variable)) +
        geom_point() +
        geom_hline(yintercept = 1) +
        geom_vline(xintercept = 7500000, color = "red") +
        theme_classic()
hits <- rownames(geneinfo)[geneinfo$start < 7500000 & geneinfo$chr == "chr17"]

# Chr22 ####
geneinfo <- (infer@gene_order)
ok <- which(geneinfo$chr == "chr22")
aframe <- data.frame(expr[ok,])
means <- do.call(cbind, lapply(infer@observation_grouped_cell_indices, function(x) rowMeans(expr[ok, x])))
means <- data.frame(means)
means$position <- geneinfo$start[ok]
tmp <- reshape2::melt(means, id.vars = "position")
ggplot(tmp, aes(position, value, group = variable, color = variable)) +
        geom_point() +
        geom_hline(yintercept = 1) +
        geom_vline(xintercept = c(18000000,30000000), color = "red") +
        theme_classic()
hits <- rownames(geneinfo)[geneinfo$start < 30000000 & geneinfo$start > 18000000 & geneinfo$chr == "chr22"]
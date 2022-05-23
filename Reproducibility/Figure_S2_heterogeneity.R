setwd("/Users/yanfangfang/Downloads/MW/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(pdist)

# Figure S2A ####
seu <- readRDS("data/MCL.rds")
meta <- seu@meta.data
ok <- which(meta$source == "PB" & meta$nCount_RNA > 1000)
seu <- seu[, ok]
meta <- meta[ok, ]
# Merge cell types into lineages ####
meta$cellgroup <- meta$celltype
meta$cellgroup[which(meta$celltype %in% c("CD8 effector", "CD8_EarlyActiv", "CD8_EffectorMemory", "CD8_NaiveLike", "CD8_Tex", "CD8_Tpex"))] <- "CD8 T"
meta$cellgroup[which(meta$celltype %in% c("NK", "NK dim", "NKT"))] <- "NK"
meta$cellgroup[which(meta$celltype %in% c("Double negative T cell", "Tfh", "Th1", "Treg", "CD4_NaiveLike"))] <- "CD4 T"
groups <- paste(meta$sample, meta$cellgroup, sep = "|")
ok <- names(which(table(groups) > 100))
tmp <- do.call(rbind, lapply(ok, function(x) strsplit(x, "|", fixed =T)[[1]]))
tmp <- data.frame(tmp)
colnames(tmp) <- c("sample", "celltype")
tmp$num_cells <- as.numeric(unclass(table(groups)[which(table(groups) > 100)]))

# Define Down-sampling function ####
Down_Sample_Matrix <- function (expr_mat)
{
        min_lib_size <- min(colSums(expr_mat))
        down_sample <- function(x) {
                prob <- min_lib_size/sum(x)
                return(unlist(lapply(x, function(y) {
                        rbinom(1, y, prob)
                })))
        }
        down_sampled_mat <- apply(expr_mat, 2, down_sample)
        return(down_sampled_mat)
}

# Calculate variance explained by patient variable ####
tmp <- tmp[tmp$sample != "Normal_1",]
tmp <- tmp[tmp$sample != "Normal_2",]
good_celltypes <- names(which(table(tmp$celltype) > 5))
calc_per_explained <- function(celltype){
        samples <- tmp$sample[tmp$celltype == celltype]
        ok <- which(meta$sample %in% samples & meta$cellgroup == celltype)
        subset <- seu[,ok]
        asplit <- split(1:ncol(subset), subset@meta.data$sample)
        cells <- unlist(lapply(asplit, function(x) sample(x, 100)))
        
        subset2 <- subset[, cells]
        counts <- subset2@assays$RNA@counts
        counts <- counts[which(apply(counts, 1, function(x) sum(x > 0)) > ncol(counts) * 0.1),]
        
        counts_ds <- Down_Sample_Matrix(counts)
        log_counts <- log(counts_ds + 1)
        
        perc_patient <- apply(log_counts, 1, function(x){
                aframe <- data.frame(gene = x, meta[colnames(log_counts),])
                afit <- try(summary(aov(gene ~ chemistry + patient, data = aframe))[[1]])
                if(class(afit) == 'try-error') return(NA)
                afit[2,2]/sum(afit[,2])
        })
}
res <- lapply(good_celltypes, calc_per_explained)
names(res) <- good_celltypes
res <- lapply(res, na.omit)
aframe <- do.call(rbind, lapply(names(res), function(x) data.frame(celltype = x, var_expl = res[[x]])))
aframe$tumor <- "no"
aframe$tumor[aframe$celltype == "Tumor B"] <- "yes"
library(ggpubr)
my_comparisons <- list(c("Tumor B","NK"),
                       c("Tumor B","CD8 T"),
                       c("Tumor B","CD4 T"),
                       c("Tumor B","CD14+ Monocytes"))
ggplot(aframe[aframe$var_expl<=0.04,], aes(celltype, var_expl, color = tumor)) +
        geom_boxplot(outlier.shape = NA) +
        #scale_y_continuous(limits = quantile(aframe$var_expl, c(0.1, 0.9))) +
        #coord_cartesian(ylim = c(0, 0.06))+
        #scale_color_manual(values = c("black", "red")) +
        ylab("Sum of squares patient/Total sum of squares") +
        labs(title = "Patient heterogeneity across cell types") +
        theme_bw() +xlab("Cell type")+
        stat_compare_means(comparisons = my_comparisons) +
        stat_compare_means(label.y = 0.04)

# Figure S2B ####
seu <- readRDS("data/3prime_w_celltype.rds")
seu <- seu[,Idents(seu)!="Eryth"]
runSeurat <- function(seu){
        seu <- subset(seu, subset = nFeature_RNA > 200 & percent.mt < 20 & nFeature_RNA < 5000)
        seu <- NormalizeData(seu)
        seu <- FindVariableFeatures(seu)
        seu <- ScaleData(seu, features = VariableFeatures(seu))
        seu <- RunPCA(seu,features=VariableFeatures(seu))
        seu <- FindNeighbors(seu, dims = 1:10)
        seu <- FindClusters(seu, resolution = 0.5)
        seu <- RunUMAP(seu, dims = 1:10)
        seu
}
tumor <- seu[,seu$celltype=="Tumor B" & seu$patient !="Normal_1" & seu$patient !="Normal_2"]
tumor <- runSeurat(tumor)
p1 <- DimPlot(tumor,group.by = 'patient',label=T,repel = T)+
        labs(title="Cohort1, Tumor B cells",subtitle="colored by patient")+NoLegend()
non_tumor <- seu[,seu$celltype!="Tumor B"]
non_tumor <- runSeurat(non_tumor)
p2 <- DimPlot(non_tumor,group.by = 'patient')+
        labs(title="Cohort1, TME cells",subtitle="colored by patient")+NoLegend()
p3 <- DimPlot(non_tumor,group.by = 'celltype')+
        labs(title="",subtitle="colored by cell type")
p1|p2|p3

seu <- readRDS("data/5prime_w_celltype.rds")
seu <- seu[,Idents(seu)!="Eryth"]
tumor <- seu[,seu$celltype=="Tumor B" & seu$patient !="Normal_1" & seu$patient !="Normal_2"]
tumor <- runSeurat(tumor)
p4 <- DimPlot(tumor,group.by = 'patient',label=T,repel = T)+
        labs(title="Cohort2, Tumor B cells",subtitle="colored by patient")+NoLegend()
non_tumor <- seu[,seu$celltype!="Tumor B"]
non_tumor <- runSeurat(non_tumor)
p5 <- DimPlot(non_tumor,group.by = 'patient')+
        labs(title="Cohort2, TME cells",subtitle="colored by patient")+NoLegend()
p6 <- DimPlot(non_tumor,group.by = 'celltype')+
        labs(title="",subtitle="colored by cell type")
p4|p5|p6

(p1|p2|p3)/(p4|p5|p6)

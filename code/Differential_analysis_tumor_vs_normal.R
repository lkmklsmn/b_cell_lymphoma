# Load R libs ####
library(Seurat)
library(ggplot2)
library(gridExtra)
library(sctransform)
library(MAST)
library(scater)

# Load data ####
setwd("/Users/lukas.simon/OneDrive/Miko/UTHealth/projects/Wang_collab/")
seu <- readRDS("data/bbknn_all_genes.rds")
meta <- seu@meta.data
meta$Therapeutic.sensitivity <- gsub("Ibrutinib naïve", "Ibrutinib naive", meta$Therapeutic.sensitivity)
meta$tumor <- "yes"
meta$tumor[grep("Normal", meta$Sample)] <- "no"
counts <- seu@assays$RNA@counts

# Restrict to PB
ok <- which(meta$source_chemistry %in% c("PB_3primer", "PB_5primer"))
meta <- meta[ok,]
counts <- counts[,ok]

# Restrict to B cells
ok <- which(meta$celltype == "B")
meta <- meta[ok,]
counts <- counts[,ok]

# Calculate pseudobulks ####
asplit <- split(rownames(meta), meta$Sample)
good <- names(which(table(meta$Sample) > 20))
asplit <- asplit[good]
means <- do.call(cbind, lapply(asplit, function(x) Matrix::rowMeans(counts[,x])))
meta_sample <- meta[match(colnames(means), meta$Sample), 1:11]
meta_sample$tumor <- "yes"
meta_sample$tumor[grep("Normal", meta_sample$Sample)] <- "no"
meta_sample$median_total_counts <- unlist(lapply(asplit, function(x) median(Matrix::colSums(counts[,x]))))

# Run pseudobulk PCA ####
ok <- which(apply(means, 1, function(x) sum(x > 1e-2)) > 10)
pca <- prcomp(t(sqrt(means[ok,])), scale. = T)
aframe <- data.frame(meta_sample, pca$x)

p_tumor <- ggplot(aframe, aes(PC1, PC2, color = tumor)) +
  geom_point() + theme_bw() + ggtitle("Pseudobulk PCA")

p_totalcounts <- ggplot(aframe, aes(PC1, PC2, color = log(median_total_counts))) +
  geom_point() + theme_bw() + ggtitle("Pseudobulk PCA")

p_chrm <- ggplot(aframe, aes(PC1, PC2, color = chemistry)) +
  geom_point() + theme_bw() + ggtitle("Pseudobulk PCA")

p_ibu <- ggplot(aframe, aes(PC1, PC2, color = Ibrutinib.sensitivity)) +
  geom_point() + theme_bw() + ggtitle("Pseudobulk PCA")

p_therapy <- ggplot(aframe, aes(PC1, PC2, color = Therapeutic.sensitivity)) +
  geom_point() + theme_bw() + ggtitle("Pseudobulk PCA")

p_cart <- ggplot(aframe, aes(PC1, PC2, color = CAR.T.senstivity)) +
  geom_point() + theme_bw() + ggtitle("Pseudobulk PCA")

grid.arrange(p_chrm, p_totalcounts, p_tumor,
             p_ibu, p_cart, p_therapy, ncol = 3)

# Find genes associated w cancer ####
three_prime <- which(meta_sample$chemistry == "3primer")
five_prime <- which(meta_sample$chemistry == "5primer")
res <- t(apply(sqrt(means[ok, three_prime]), 1, function(x){
  coefficients(summary(lm(x ~ meta_sample$tumor[three_prime])))[2, c(1, 4)]
}))
res <- data.frame(res)
colnames(res) <- c("estimate", "pval")
res$label <- NA
cancer_genes <- c("CCND1", "MYC", "SOX11", "BTK")
res[cancer_genes, "label"] <- cancer_genes

ggplot(res, aes(estimate, -log10(pval), label = label)) + 
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = 'log', breaks = c(1, 10, 100)) +
  geom_label() +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-2, 2) +
  ylab("P-value [-log10]") + xlab("Estimate") + ggtitle("Pseudobulk: cancer vs normal (3 primer only)") +
  theme_bw()

good_genes <- rownames(res)[which(abs(res$estimate) > 0.2 & res$pval < 0.1)]

# Create Seurat object ####
ok <- which(meta$chemistry == "3primer")
subset <- CreateSeuratObject(counts[,ok])
subset@meta.data <- meta[ok,]
subset <- SCTransform(subset, verbose = 0)
subset <- ScaleData(subset, features = good_genes, verbose = 0)
subset <- RunPCA(subset, features = good_genes, verbose = 0)
plot(subset@reductions$pca@stdev, ylab = "Std Dev", xlab = "PCs")

subset <- FindNeighbors(subset, dims = 10, verbose = 0)
subset <- RunUMAP(subset, dims = 1:10, verbose = 0)
DimPlot(object = subset, reduction = 'umap', group.by = "tumor")
DimPlot(object = subset, reduction = 'pca', group.by = "tumor")
DimPlot(object = subset, reduction = 'umap', group.by = "Sample", label = T)
DimPlot(object = subset, reduction = 'umap', group.by = "patient")
DimPlot(object = subset, reduction = 'umap', group.by = "Ibrutinib.sensitivity")

FeaturePlot(object = subset, features = cancer_genes, reduction = 'umap')

# Run cancer vs normal DE analysis ####
sce <- as.SingleCellExperiment(subset)
sca <- SceToSingleCellAssay(sce)

asplit <- split(rownames(sca@colData), sca@colData$Sample)
cells <- unlist(lapply(asplit, function(x) sample(x, min(length(x), 20))))

tmp <- sca[,cells]
ok <- names(which(apply(assays(tmp)[["counts"]], 1, function(x) sum(x > 0)) > 50))

res <- zlm( ~ tumor + (1|tumor:patient),
            sca = tmp[ok,],
            exprs_value = 'logcounts',
            method='glmer', ebayes=FALSE)
tmp <- summary(res)$datatable


# 

library(limma)

genes <- names(which(rowMeans(counts[, cells]) > 0.1))

design <- model.matrix( ~ tumor, sca@colData[cells,])
vobj <- voom(counts[genes, cells], design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj, design, block = sca@colData[cells, "patient"])
fitDupCor <- lmFit(vobj, design, block = sca@colData[cells, "patient"], correlation = dupcor$consensus)
fitDupCor <- eBayes(fitDupCor)
res <- topTable(fitDupCor, number = nrow(vobj))

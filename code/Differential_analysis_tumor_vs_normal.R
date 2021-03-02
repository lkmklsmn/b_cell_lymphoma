# Load R libs ####
library(Seurat)
library(ggplot2)
library(gridExtra)
library(sctransform)
library(MAST)
library(scater)
library(pheatmap)

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

# Run cancer vs normal DE analysis using MAST ####
sce <- as.SingleCellExperiment(subset)
sca <- SceToSingleCellAssay(sce)

# Subsample to a maximum of 100 cells per patient
asplit <- split(rownames(sca@colData), sca@colData$Sample)
good <- names(which(unlist(lapply(asplit, length)) > 20))
asplit <- asplit[good]
cells <- unlist(lapply(asplit, function(x) sample(x, min(length(x), 100))))

# Subset to genes detected in more than 10% of cells in at least 2 samples
freq_detected <- do.call(cbind, lapply(asplit, function(x) apply(assays(sca)[["counts"]][, x], 1, function(y) mean(y > 0))))
ok <- which(apply(freq_detected, 1, function(x) sum(x > 0.2)) >= 2)

# Run mixed model with random effect (taken from https://github.com/kdzimm/PseudoreplicationPaper/blob/c3059a3b361e89bde595f222757d04b89f77eb62/Power/Power.Rmd#L60)
# takes ~ 40min
system.time(res <- suppressMessages(
  zlm( ~ tumor + (1|patient),
       sca = sca[ok, cells],
       exprs_value = 'logcounts',
       strictConvergence = F, method='glmer', ebayes=F)))

# Format results into single table
# Takes ~25min
system.time(summaryCond <- suppressMessages(summary(res, doLRT='tumoryes')))
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='tumoryes'
                            & summaryDt$component=='logFC', c(1,7,5,6,8)],
                  summaryDt[summaryDt$contrast=='tumoryes'
                            & summaryDt$component=='H', c(1,4)],
                  by = 'primerid')
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle <- fcHurdle[sort.list(fcHurdle$`Pr(>Chisq)`),]

# Heatmap of significant genes ####
sig <- fcHurdle$primerid[which(fcHurdle$`Pr(>Chisq)` < 1e-3 & abs(fcHurdle$z) > 4)]
tmp <- data.matrix(subset@assays$SCT@data[sig, cells])
tmp <- tmp[, sample(colnames(tmp), 200)]
pheatmap(tmp, annotation_col = meta[, c("tumor", "patient")],
         show_rownames = F, show_colnames = F,
         scale = "row", breaks = seq(-2, 2, length = 101))

ggplot(fcHurdle, aes(x = z, y = -log10(`Pr(>Chisq)`))) +
  geom_point() + theme_bw()

flat_dat <- as(sca[sig[1:6],], 'data.table')
ggplot(flat_dat, aes(x=tumor, y=logcounts,color=tumor)) +
  facet_wrap(~primerid, scale='free_y') + 
  ggtitle("DE Genes tumor vs normal") +
  geom_violin() + theme_bw()
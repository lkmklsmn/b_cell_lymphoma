# Load R libs ####
library(Seurat)
library(pheatmap)
library(ggplot2)
library(sctransform)

# Load data ####
seu <- readRDS("/Users/lukas/OneDrive/Miko/THINC/projects/Wang_collab/b_seu.rds")
res <- read.delim("/Users/lukas/OneDrive/Miko/THINC/projects/Wang_collab/mixed_model_results.txt")
specific_genes <- lapply(c("dual", "r", "s", "slow"), function(x) tail(rownames(res[sort.list(res[, paste(x, "spe", sep = "_")]),]), 200))

# Create pseudobulk heatmap ####
asplit <- split(rownames(seu@meta.data), seu@meta.data$sample)
ok <- names(which(unlist(lapply(asplit, length)) > 10))
asplit <- asplit[ok]
means <- do.call(cbind, lapply(asplit, function(x) Matrix::rowMeans(seu@assays$RNA@data[, x])))
colnames(means) <- names(asplit)

anno <- seu@meta.data[match(colnames(means), seu@meta.data$sample), c("ibrutinib_sensitivity", "patient")]
rownames(anno) <- colnames(means)

pheatmap(means[unlist(specific_genes), order(anno$ibrutinib_sensitivity, anno$patient)],
         scale = "row",
         cluster_rows = F, cluster_cols = F,
         annotation_col = anno,
         show_rownames = F, show_colnames = F,
         breaks = seq(-2, 2, length = 101))

# Create PCA plots ####
sub <-subset(seu, cells = rownames(seu@meta.data[seu@meta.data$sample %in% ok,]))

sub <- ScaleData(sub, assay = "RNA", features = specific_genes[[1]])
sub <- RunPCA(sub, features = specific_genes[[1]], assay = "RNA", npcs = 4)
aframe <- data.frame(sub@meta.data, sub@reductions$pca@cell.embeddings)
aframe$sample <- factor(aframe$sample, levels = unique(aframe$sample[order(aframe$ibrutinib_sensitivity)]))
p_dual <- ggplot(aframe, aes(sample, PC_1, color = ibrutinib_sensitivity)) +
  geom_boxplot() +
  ggtitle("Dual-specific genes") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
sub <- ScaleData(sub, assay = "RNA", features = specific_genes[[2]])
sub <- RunPCA(sub, features = specific_genes[[2]], assay = "RNA", npcs = 4)
aframe <- data.frame(sub@meta.data, sub@reductions$pca@cell.embeddings)
aframe$sample <- factor(aframe$sample, levels = unique(aframe$sample[order(aframe$ibrutinib_sensitivity)]))
p_r <- ggplot(aframe, aes(sample, PC_1, color = ibrutinib_sensitivity)) +
  geom_boxplot() +
  ggtitle("R-specific genes") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sub <- ScaleData(sub, assay = "RNA", features = specific_genes[[3]])
sub <- RunPCA(sub, features = specific_genes[[3]], assay = "RNA", npcs = 4)
aframe <- data.frame(sub@meta.data, sub@reductions$pca@cell.embeddings)
aframe$sample <- factor(aframe$sample, levels = unique(aframe$sample[order(aframe$ibrutinib_sensitivity)]))
p_s <- ggplot(aframe, aes(sample, PC_1, color = ibrutinib_sensitivity)) +
  geom_boxplot() +
  ggtitle("S-specific genes") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

sub <- ScaleData(sub, assay = "RNA", features = specific_genes[[4]])
sub <- RunPCA(sub, features = specific_genes[[4]], assay = "RNA", npcs = 4)
aframe <- data.frame(sub@meta.data, sub@reductions$pca@cell.embeddings)
aframe$sample <- factor(aframe$sample, levels = unique(aframe$sample[order(aframe$ibrutinib_sensitivity)]))
p_slow <- ggplot(aframe, aes(sample, PC_1, color = ibrutinib_sensitivity)) +
  geom_boxplot() +
  ggtitle("Slow responder-specific genes") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

gridExtra::grid.arrange(p_dual, p_r, p_s, p_slow)

# #####
combined <- readRDS("/Users/lukas/OneDrive/Miko/THINC/projects/Wang_collab/combined.rds")
combined_b_cells <- subset(combined,idents = c(6,16,10,11,1,4,12))
fiveprime <- subset(combined_b_cells, cells = rownames(combined@meta.data[combined@meta.data$CART_sensitivity %in% c("R", "S"), ]))

sub <- ScaleData(fiveprime, features = unlist(specific_genes))
sub <- RunPCA(sub, features = specific_genes[[1]], npcs = 10)
sub <- RunUMAP(sub, dims = 1:10)

sub <- SCTransform(fiveprime, vars.to.regress = "percent.mt", verbose = FALSE)
sub <- ScaleData(sub, assay = "SCT", features = unlist(specific_genes))  
  
sub@meta.data$dual_sig <- colMeans(sub@assays$SCT@scale.data[specific_genes[[1]],])
sub@meta.data$r_sig <- colMeans(sub@assays$SCT@scale.data[specific_genes[[2]],])

DimPlot(sub, reduction = "umap", group.by = "CART_sensitivity")
FeaturePlot(sub, features = c("dual_sig", "r_sig"), reduction = "umap")
VlnPlot(sub, features = c("dual_sig", "r_sig"), group.by = "CART_sensitivity")
VlnPlot(sub, features = c("nCount_RNA"), group.by = "CART_sensitivity")

ok <- intersect(rownames(seu), unlist(specific_genes[1:2]))
ok <- intersect(ok, rownames(sub))
seu <- ScaleData(seu, features = ok)
cells <- rownames(seu@meta.data[seu@meta.data$ibrutinib_sensitivity %in% c("Dual", "R"),])
pca <- prcomp(t(seu@assays$RNA@scale.data[ok, cells]))

aframe <- data.frame(seu@meta.data[cells, ], pca$x)
ggplot(aframe, aes(PC1, PC2, color = ibrutinib_sensitivity)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  theme_bw()

sub <- ScaleData(sub, assay = "RNA", features = ok)
pca2 <- predict(pca, t(sub@assays$RNA@scale.data[ok, ]))

aframe <- data.frame(sub@meta.data, pca2)
ggplot(aframe, aes(PC1, PC2, color = CART_sensitivity)) +
  geom_point() +
  xlab("PC1") + ylab("PC2") +
  theme_bw()

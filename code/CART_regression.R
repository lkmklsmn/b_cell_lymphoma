args<-commandArgs(T)

# Load R libs ####
library(Seurat)
library(ggplot2)
library(gridExtra)
library(sctransform)
library(MAST)
library(scater)
library(pheatmap)
library(reshape2)

# Load data ####
seu <- readRDS(args[1])
meta <- seu@meta.data
counts <- seu@assays$RNA@counts

print("Restrict to CART samples")
ok <- which(meta$chemistry == "5primer")
meta <- meta[ok,]
counts <- counts[,ok]

print("Restrict to B cells")
ok <- which(meta$celltype == "B")
meta <- meta[ok,]
counts <- counts[,ok]

print("Create Seurat object")
subset <- CreateSeuratObject(counts)
subset@meta.data <- meta
subset <- SCTransform(subset, verbose = 0)
subset <- ScaleData(subset, verbose = 0)

print("Run CART R vs S DE analysis using MAST")
sce <- as.SingleCellExperiment(subset)
sca <- SceToSingleCellAssay(sce)

# Subsample to a maximum of 200 cells per sample
asplit <- split(rownames(sca@colData), sca@colData$Sample)
good <- names(which(unlist(lapply(asplit, length)) > 20))
asplit <- asplit[good]
cells <- unlist(lapply(asplit, function(x) sample(x, min(length(x), 200))))

# Subset to genes detected in more than 10% of cells in at least 2 samples
freq_detected <- do.call(cbind, lapply(asplit, function(x) apply(assays(sca)[["counts"]][, x], 1, function(y) mean(y > 0))))
ok <- names(which(apply(freq_detected, 1, function(x) sum(x > 0.2)) >= 2))

print("Run mixed model with random effect")
system.time(res <- suppressMessages(
  zlm( ~ CAR.T.senstivity + (1|Sample),
       sca = sca[ok, cells],
       exprs_value = 'logcounts',
       strictConvergence = F, method='glmer', ebayes=F)))

print("Format results into single table")
system.time(summaryCond <- suppressMessages(summary(res, doLRT='CAR.T.senstivityS')))
summaryDt <- summaryCond$datatable
fcHurdle <- merge(summaryDt[summaryDt$contrast=='CAR.T.senstivityS'
                            & summaryDt$component=='logFC', c(1,7,5,6,8)],
                  summaryDt[summaryDt$contrast=='CAR.T.senstivityS'
                            & summaryDt$component=='H', c(1,4)],
                  by = 'primerid')
fcHurdle <- stats::na.omit(as.data.frame(fcHurdle))
fcHurdle$fdr <- p.adjust(fcHurdle$`Pr(>Chisq)`, 'fdr')
fcHurdle <- fcHurdle[sort.list(fcHurdle$`Pr(>Chisq)`),]
saveRDS(fcHurdle,file="/home/fyan/CART_fcHurdle.rds")

print("print plots")
load("/home/fyan/BlueYellowColormaps_V1.RData")
sig <- fcHurdle$primerid[which(fcHurdle$`Pr(>Chisq)` < 1e-3 & abs(fcHurdle$z) > 2)]
reps <- 3
tmp <- data.matrix(subset@assays$SCT@data[sig, ])
means <- do.call(cbind, lapply(asplit, function(x){
  cells <- lapply(1:3, function(y) sample(x, 10))
  do.call(cbind, lapply(cells, function(y) rowMeans(tmp[,y])))
}))
colnames(means) <- unlist(lapply(names(asplit), function(x) paste(x, 1:3, sep = "_")))
anno <- data.frame(response = rep("CART_R", ncol(means)))
anno$sample <- substr(colnames(means),1,2)
anno$response <- subset@meta.data[match(anno$sample,subset@meta.data$Sample),
                                  'CAR.T.senstivity']
rownames(anno) <- colnames(means)
pdf(file = "/home/fyan/CART_heatmap.pdf",width = 6,height = 4)
pheatmap(means, annotation_col = anno,
         show_rownames = F, show_colnames = T,
         scale = "row", breaks = seq(-2, 2, length = 257),
         color = yellow2blue, main = 'CAR.T R vs S')
dev.off()

# Volcano plot of differentially expressed genes ####
pdf(file = "/home/fyan/CART_volcano.pdf",width = 6,height = 4)
ggplot(fcHurdle, aes(x = coef, y = -log10(`Pr(>Chisq)`))) +
  geom_hex(bins = 50) +
  scale_fill_viridis_c(trans = 'log', breaks = c(1, 10, 100)) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  xlim(-1, 1) + ggtitle("Volcano: ~ CAR.T.senstivity + (1|Sample") +
  theme_bw()
dev.off()

# Violin plots ####
sig <- fcHurdle$primerid[1:12]
flat_dat <- data.frame(CART = "R", 
                       sample = colnames(freq_detected),
                       t(freq_detected[sig,]))
flat_dat$CART <- subset@meta.data[match(flat_dat$sample,subset@meta.data$Sample),
                 'CAR.T.senstivity']
flat_dat <- melt(flat_dat, id.vars = c("CART", "sample"))
pdf(file = "/home/fyan/CART_violin.pdf",width = 6,height = 4)
ggplot(flat_dat, aes(x=CART, y=value,color=CART)) +
  facet_wrap(~variable, scale='free_y') + 
  ggtitle("DE Genes, CAR.T R vs S") +
  geom_boxplot() + geom_point() +
  theme_bw()
dev.off()

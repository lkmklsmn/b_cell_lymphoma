# Load R libs ####
library(Seurat)
library(nlme)

# Load data ####
seu <- readRDS("/Users/lukas/Downloads/b_seu.rds")

# Plot tSNE embedding ####
DimPlot(seu, group.by = "sample", label = T)
DimPlot(seu, group.by = "patient", label = T)

# Define expression and sample annotation objects ####
meta <- seu@meta.data
meta$ibrutinib_sensitivity <- factor(meta$ibrutinib_sensitivity, levels = c("Normal", "Dual", "R", "S", "Slow_responder"))
meta$patient <- factor(meta$patient, levels = c("Normal", "A", "L", "E", "W", "Y", "B", "C", "D", "V", "AA", "X", "Z"))
meta$sample <- factor(meta$sample, levels = unique(meta$sample[order(meta$patient)]))
expr <- seu@assays$RNA@data

# Define genes for analysis ####
sums <- log(rowSums(seu@assays$RNA@counts))
genes <- names(which(sums > 3))

# Run regression for 12k genes & all cells takes about 35min ####
system.time({
  res <- do.call(rbind, lapply(genes, function(x){
    #print(x)
    aframe <- data.frame(meta, expr = expr[x,])
    #aframe <- aframe[r,]
    m1 <- try(lme(expr ~ S.Score + G2M.Score + ibrutinib_sensitivity, random = ~ 1|patient, data = aframe))
    if(class(m1) == "try-error") return(rep(NA, 8))
    as.numeric(coefficients(summary(m1))[4:7, c(1, 5)])
  }))
})
rownames(res) <- genes
colnames(res) <- c("coef_Dual", "coef_R", "coef_S", "coef_slow", "pval_Dual", "pval_R", "pval_S", "pval_slow")
res <- data.frame(res)

#Define plotting function ####
create_plot <- function(gene = "CD52"){
  aframe <- data.frame(meta, expr = expr[gene,])
  ggplot(aframe, aes(patient, expr, group = sample, color = ibrutinib_sensitivity)) +
    geom_boxplot() +
    ggtitle(gene) + ylab("Expression levels") +
    theme_bw()
}
create_plot()

# Show some examples ####
FeaturePlot(seu, features = "PRPF19")
create_plot("PRPF19")

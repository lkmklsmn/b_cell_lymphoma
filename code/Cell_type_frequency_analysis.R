# Load R libs ####
library(DirichletReg)
library(ggplot2)
library(lme4)
library(reshape2)

# Load data ####
aframe <- read.delim("/Users/lukas/Downloads/cluster.txt")

# Extract counts and proportions ####
cell_type_counts <- unclass(table(aframe$sample, aframe$functional.cluster))
meta <- aframe[match(rownames(cell_type_counts), aframe$sample), c("sample", "patient", "source", "ibrutinib_sensitivity", "chemistry")]
cell_type_prop <- cell_type_counts/rowSums(cell_type_counts)

# Plot some proportions ####
tmp <- data.frame(cell_type_prop, meta)
tmp <- melt(tmp, measure.vars = colnames(cell_type_prop))
ggplot(tmp, aes(patient, sqrt(value), color = chemistry)) +
  facet_wrap(~ variable) +
  geom_boxplot() + geom_point() +
  theme_bw()


# Plot total number of T cells ####
ggplot(aframe, aes(ibrutinib_sensitivity, total_counts, color = ibrutinib_sensitivity)) +
  facet_wrap(~ chemistry, scales = "free_x") +
  geom_boxplot() +
  theme_bw()

# Dimension reduced view ####
ok <- which(meta$chemistry == "3prime")
d <- dist(sqrt(cell_type_prop[ok, ]))
fit <- cmdscale(d,eig=TRUE, k=2)
tmp <- data.frame(fit$points, meta[ok, ])

p3 <- ggplot(tmp, aes(X1, X2, color = ibrutinib_sensitivity )) +
  geom_point() +
  stat_ellipse() +
  ggtitle("MDS of T cell proportions - 3prime") +
  theme_bw()

ok <- which(meta$chemistry == "5prime")
d <- dist(sqrt(cell_type_prop[ok, ]))
fit <- cmdscale(d,eig=TRUE, k=2)
tmp <- data.frame(fit$points, meta[ok, ])

p5 <- ggplot(tmp, aes(X1, X2, color = ibrutinib_sensitivity )) +
  geom_point() +
  stat_ellipse() +
  ggtitle("MDS of T cell proportions - 5prime") +
  theme_bw()

gridExtra::grid.arrange(p3, p5, ncol = 2)

# Regression ####
# Binomial approach 
total_counts <- rowSums(cell_type_counts)
pvals_table <- t(apply(cell_type_counts, 2, function(x){
  counts <- cbind(x, total_counts)
  #afit <- glm(counts ~ meta$chemistry + meta$ibrutinib_sensitivity, family = binomial(link = "logit"))
  #coefficients(summary(afit))[-c(1:2), 4]
  afit <- glmer(counts ~ meta$chemistry + (1 | meta$patient) + meta$ibrutinib_sensitivity, family = binomial(link = "logit"))
  coefficients(summary(afit))[-c(1:2), 4]
}))
tmp <- -log10(pvals_table)
colnames(tmp) <- gsub("meta$ibrutinib_sensitivity", "", colnames(tmp), fixed = T)
pheatmap::pheatmap(tmp, cluster_rows = F, cluster_cols = F,
                   main = "Significant associations binomial model")

# Define plotting function ####
plot_celltype <- function(cell_type="CD8_EffectorMemory"){
  tmp <- data.frame(proportion = cell_type_prop[, cell_type], meta)
  ggplot(tmp, aes(ibrutinib_sensitivity, proportion, color = ibrutinib_sensitivity)) +
    facet_wrap(~ chemistry, scales = "free_x") +
    geom_boxplot() + geom_point() +
    ylab(paste("Proportion", cell_type)) +
    theme_bw()
}


# Dirichlet approach ####
aframe <- meta
aframe$ibrutinib_sensitivity <- factor(aframe$ibrutinib_sensitivity, levels = c("Normal", "S", "Slow_responder", "R", "Dual"))
ok <- which(aframe$chemistry == "3prime")
aframe <- aframe[ok,]
aframe$Y <- DR_data(cell_type_counts[ok,])
afit <- DirichReg(Y ~ source + ibrutinib_sensitivity, aframe)
summary(afit)

aframe <- meta
aframe$ibrutinib_sensitivity <- factor(aframe$ibrutinib_sensitivity, levels = c("Normal", "S", "Slow_responder", "R", "Dual"))
ok <- which(aframe$chemistry == "5prime")
aframe <- aframe[ok,]
aframe$Y <- DR_data(cell_type_counts[ok,])
afit <- DirichReg(Y ~ ibrutinib_sensitivity, aframe)
summary(afit)

gridExtra::grid.arrange(plot_celltype("CD8_EffectorMemory"),
                        plot_celltype("Th1"), ncol = 2)

aframe <- meta
aframe$ibrutinib_sensitivity <- factor(aframe$ibrutinib_sensitivity, levels = c("Normal", "S", "Slow_responder", "R", "Dual"))
aframe$Y <- DR_data(cell_type_counts)
afit <- DirichReg(Y ~ chemistry + ibrutinib_sensitivity, aframe)
summary(afit)


# Create plots ####
# Best in 5prime (dual vs R)
p1 <- plot_celltype("Th1")
p2 <- plot_celltype("CD8_Tex")
p3 <- plot_celltype("Treg")
gridExtra::grid.arrange(p1, p2, p3, ncol = 3)

plot_celltype("CD8_EffectorMemory")

# Best in 3prime 
plot_celltype("CD8_EffectorMemory")
plot_celltype("Treg")
plot_celltype("CD8_Tex")
plot_celltype("CD8_NaiveLike")
plot_celltype("CD8_EarlyActiv")

# Load R libs 
library(MAST)
library(scater)
library(Seurat)
library(sctransform)

# Load data
seu <- readRDS("/Users/lukas.simon/OneDrive/Miko/UTHealth/projects/Wang_collab/data/bbknn.rds")

meta <- seu@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1", "Normal_2"))] <- "no"

counts <- seu@assays$RNA@counts

ok <- which(meta$celltype == "B")
meta <- meta[ok,]
counts <- counts[,ok]

ok <- which(meta$Sample.source == "PB")
meta <- meta[ok,]
counts <- counts[,ok]

seu2 <- CreateSeuratObject(counts = counts)
seu2@meta.data <- meta
seu2 <- SCTransform(seu2, vars.to.regress = "chemistry")


# Load data
seu <- readRDS("/Users/lukas.simon/OneDrive/Miko/UTHealth/projects/Wang_collab/data/bbknn.rds")

meta <- seu@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1", "Normal_2"))] <- "no"

counts <- seu@assays$RNA@counts

ok <- which(meta$celltype == "B")
meta <- meta[ok,]
counts <- counts[,ok]

ok <- which(meta$Sample.source == "PB")
meta <- meta[ok,]
counts <- counts[,ok]

asplit <- split(rownames(meta), meta$patient)
r <- unlist(lapply(asplit, function(x) sample(x, min(100, length(x)))))

counts <- counts[,r]
meta <- meta[r,]

r <- sample(rownames(counts), 5000)
r <- c("CCND1", r)
#r <- c("CCND1", "UBE2QL1")
counts <- data.matrix(counts[r,])

good <- names(which(table(meta$patient) > 26))
ok <- which(meta$patient %in% good)
meta <- meta[ok,]
meta$patient <- as.factor(as.character(meta$patient))
counts <- counts[,ok]

sce <- SingleCellExperiment(
  assays = list(counts = counts), 
  colData = meta
)
sce

sce <- logNormCounts(sce)
sce <- addPerCellQC(sce)
plotColData(sce, y="detected", x="total")

sce <- runPCA(sce, ncomponents=5, exprs_values = 'logcounts')
plotReducedDim(sce, dimred = 'PCA', colour_by = 'chemistry')

sca <- SceToSingleCellAssay(sce)
res <- zlm( ~ tumor + (1|tumor:patient),
            sca = sca,
            exprs_value = 'logcounts',
            method='glmer', ebayes=FALSE)
res <- summary(res)$datatable
tmp <- res
tmp <- na.omit(tmp)

tmp <- tmp[tmp$component == "D" & tmp$contrast == "tumoryes",]
tmp <- tmp[rev(sort.list(abs(tmp$z))),]
head(tmp)

plot_gene <- function(gene = "OAZ1"){
  expr_gene <- assays(sca)[["logcounts"]][gene,]
  aframe <- data.frame(expr_gene, sca@colData)
  
  ggplot(aframe, aes(tumor, expr_gene, color = patient)) +
    geom_boxplot() + #geom_point() +
    theme_bw() + ggtitle(gene)  
}
plot_gene("OAZ1")
plot_gene("LAPTM5")
plot_gene("TMSB10")
plot_gene("RASGRP4")

# Pseudobulks ####
library(readxl)

seu <- readRDS("/Users/lukas.simon/OneDrive/Miko/UTHealth/projects/Wang_collab/data/bbknn.rds")

meta <- seu@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1", "Normal_2"))] <- "no"

counts <- seu@assays$RNA@counts

ok <- which(meta$celltype == "B")
meta <- meta[ok,]
counts <- counts[,ok]

ok <- which(meta$Sample.source == "PB")
meta <- meta[ok,]
counts <- counts[,ok]

asplit <- split(rownames(meta), meta$Sample)
good <- names(which(table(meta$Sample) > 10))
asplit <- asplit[good]
ps_bulk <- do.call(cbind, lapply(asplit, function(x) rowSums(counts[,x])))

sample_info <- read_excel("/Users/lukas.simon/OneDrive/Miko/UTHealth/projects/Wang_collab/data/meta_info.xlsx", sheet = 1)
sample_info <- data.frame(sample_info)
rownames(sample_info) <- sample_info$Sample...2
sample_info <- sample_info[colnames(ps_bulk),]
sample_info$tumor <- "yes"
sample_info$tumor[which(sample_info$Sample...2 %in% c("Normal_1", "Normal_2"))] <- "no"
sample_info$experiment <- sample_info$Sample...2

library(edgeR)
library(lme4)
library(limma)


# Get total counts ####
asplit <- split(colnames(counts), meta$Sample)
good <- names(which(table(meta$Sample) > 10))
total_counts <- unlist(lapply(asplit[good], function(x) median(Matrix::colSums(counts[,x]))))
sample_info$total_counts <- total_counts[rownames(sample_info)]

sample_info$chemistry <- meta$chemistry[match(rownames(sample_info), meta$Sample)]

gene <- "CCND1"
tmp <- table(counts[gene, ] > 0, meta$Sample)
tmp <- tmp[,rownames(sample_info)]
tmp <- t(tmp)
gene_counts <- tmp

asplit <- split(rownames(meta), meta$Sample)
good <- names(which(table(meta$Sample) > 10))
asplit <- asplit[good]
ps_bulk <- do.call(cbind, lapply(asplit, function(x) rowSums(counts[,x])))

asplit <- split(rownames(meta), meta$Sample)
good <- names(which(table(meta$Sample) > 10))
asplit <- asplit[good]
num_of_cells <- unlist(lapply(asplit, length))

asplit <- split(rownames(meta), meta$Sample)
good <- names(which(table(meta$Sample) > 10))
asplit <- asplit[good]
percent_detected <- do.call(cbind, lapply(asplit, function(x) apply(counts[,x], 1, function(y) mean(y > 0))))

ps_bulk <- ps_bulk[rowSums(ps_bulk) > 20,]

res <- do.call(rbind, lapply(rownames(ps_bulk), function(gene){
  aframe <- data.frame(expr = ps_bulk[gene,], sample_info)
  m <- try(glmer.nb(expr ~ offset(log(total_counts)) + chemistry + tumor + (1 | Patient), aframe))
  if(class(m) == "try-error") return(c(NA, NA))
  coefficients(summary(m))[3, c(1, 4)]
}))
rownames(res) <- rownames(ps_bulk)
res <- res[sort.list(res[,2]),]

per_detected_normal <- rowMeans(percent_detected[,c("Normal_1", "Normal_2")])
per_detected_cancer <- rowMeans(percent_detected[,setdiff(colnames(per_detected), c("Normal_1", "Normal_2"))])

out <- data.frame(res,
                  normal = per_detected_normal[rownames(res)],
                  cancer = per_detected_cancer[rownames(res)])
colnames(out) <- c("estimate", "pval", "normal", "cancer")
out$diff <- out$cancer - out$normal

res_lm <- do.call(rbind, lapply(rownames(ps_bulk), function(gene){
  aframe <- data.frame(expr = percent_detected[gene,], sample_info)
  afit <- lm(expr ~ num_of_cells + log(total_counts) + chemistry + tumor, data = aframe)
  coefficients(summary(afit))["tumoryes", c(1, 4)]
}))
rownames(res_lm) <- rownames(ps_bulk)
res_lm <- res_lm[sort.list(res_lm[,2]),]

head(out[which(out$estimate > 0 & out$cancer > 0.1),])




aframe <- data.frame(expr = ps_bulk[gene,], sample_info)
m <- glmer.nb(expr ~ offset(log(total_counts)) + chemistry + Ibrutinib.sensitivity + (1 | Patient),
              aframe[aframe$Ibrutinib.sensitivity %in% c("R", "S"),])

f1 <- "expr ~ offset(log(total_counts)) + chemistry + (1 | Patient)"
f2 <- "expr ~ offset(log(total_counts)) + chemistry + (1 | tumor/Patient)"
mres1 <- glmer.nb(f1, data = aframe)
mres2 <- glmer.nb(f2, data = aframe)
anova(mres1, mres2)

aframe <- data.frame(expr = counts[gene,], meta)
f1 <- "expr ~ offset(log(total_counts)) + chemistry + (1 | patient)"
f2 <- "expr ~ offset(log(total_counts)) + chemistry + (1 | tumor/patient)"
mres1 <- glmer.nb(f1, data = aframe)
mres2 <- glmer.nb(f2, data = aframe)
anova(mres1, mres2)


# Find genes up/down regulated in cancer vs normal #### 
res <- do.call(rbind, lapply(rownames(counts)[1:1000], function(gene){
  tmp <- table(counts[gene, ] > 0, meta$Sample)
  tmp <- tmp[,rownames(sample_info)]
  gene_counts <- t(tmp)
  
  m1 <- try(glmer(gene_counts ~ log(sample_info$total_counts) + sample_info$chemistry + sample_info$tumor + (1|sample_info$Patient), family = binomial(link = "logit")))
  if(class(m1) == "try-error") return(c(NA, NA))
  
  coefficients(summary(m1))[4, c(1, 4)]
}))
rownames(res) <- rownames(counts)[1:1000]
res <- res[sort.list(res[,2]),]

# Find genes up/down regulated in cancer vs rest #### 
ok <- names(which(rowSums(counts) > 10))
res <- do.call(rbind, lapply(ok, function(gene){
  tmp <- table(counts[gene, ] > 0, meta$Sample)
  tmp <- tmp[,rownames(sample_info)]
  tmp <- t(tmp)
  gene_counts <- tmp
  
  if(sum(gene_counts[,1] > 0) <= 2) return(c(NA, NA))
  
  m1 <- try(glmer(gene_counts ~ log(total_counts) + sample_info$chemistry + sample_info$tumor + (1 | sample_info$Patient), family = binomial(link = "logit")))
  if(class(m1) == "try-error") return(c(NA, NA))
  coefficients(summary(m1))[4, c(1, 4)]
}))
rownames(res) <- ok
res <- res[sort.list(res[,2]),]

gene <- "MYC"
aframe <- data.frame(expr = counts[gene,], meta)
m1 <- glmer.nb(expr ~ offset(log(total_counts)) + chemistry + tumor + (1 | patient), aframe)
m2 <- glmer.nb(expr ~ log(total_counts) + chemistry + tumor + (1 | patient), aframe)
m3 <- glmer.nb(expr ~ offset(log(total_counts)) + chemistry +  (1 + tumor | patient), aframe)

help(SCTransform)

summary(glmer(tmp ~ sample_info$Ibrutinib.sensitivity + (1|sample_info$Patient), family = binomial(link = "logit")))
summary(m2)
summary(m1)
norm <- cpm(ps_bulk)
norm <- log(norm + 1) 

v <- voom(ps_bulk)

aframe <- data.frame(expr = v$E["CCND1",], sample_info)

afit <- lm(expr ~ tumor, data = aframe)
summary(afit)

ggplot(aframe, aes(tumor, log(expr + 1), color = Sample...2)) + geom_point()

plot_gene("CCND1")

# Load R libs ####
library(Seurat)
library(sctransform)
library(lme4)
library(lmerTest)

# Load data ####
seu <- readRDS("/Users/lukas.simon/Downloads/b_cell.rds")

# Normalize data ####
seu <- SCTransform(seu)

# Define model for diff analysis ####
meta <- seu@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1", "Normal_2"))] <- "no"
meta$batch <- as.factor(as.character(meta$batch))
expr <- seu@assays$SCT@scale.data

# Visualize hits ####
gene <- "CASC7"
gene <- "CCND1"
expr_gene <- expr[gene,]

aframe <- data.frame(expr_gene, meta)
afit <- lme(expr_gene ~ tumor, random =  ~ 1|patient, data = aframe)
summary(afit)

afit1 <- lmer(expr_gene ~ patient + tumor, data = aframe)
afit2 <- lmer(expr_gene ~ tumor + (tumor|patient), data = aframe)
anove(afit1, afit2)


# Run model ####
r <- sample(rownames(expr), 500)
res <- unlist(lapply(r, function(gene){
  expr_gene <- expr[gene,]
  aframe <- data.frame(expr_gene, meta)
  #afit <- lme(expr_gene ~ tumor, random =  ~ 1|patient, data = aframe)
  afit <- lmer(expr_gene ~ tumor + (tumor|patient), data = aframe)
  summary(afit)$coefficients[2, 5]
}))
names(res) <- r

head(sort(res))

# Define plotting function ####
plot_gene <- function(gene = "CASC7"){
  expr_gene <- expr[gene,]
  aframe <- data.frame(expr_gene, meta)
  
  ggplot(aframe, aes(tumor, expr_gene, color = patient)) +
    geom_boxplot() + #geom_point() +
    theme_bw() + ggtitle(gene)  
}
plot_gene("CASC7")
plot_gene("CCND1")

#  ####
library(mgcv)
pvals_time <- do.call(rbind, lapply(r, function(gene){
  expr_gene <- expr[gene,]
  aframe <- data.frame(expr_gene, meta)
  aframe <- aframe[which(aframe$patient %in% c("V", "C")),]
  aframe$Ibrutinib.sensitivity <- as.factor(as.character(aframe$Ibrutinib.sensitivity))
  
  afit <- try(gam(expr_gene ~ Days.of.Ibrutinib.Treatment * Ibrutinib.sensitivity, data = aframe))
  
  if(class(afit) == "try-error") return(c(NA, NA, NA))
  summary(afit)[[4]][-1]
}))
rownames(pvals_time) <- r

head(pvals_time[sort.list(pvals_time[,1]), ])

head(pvals_time[sort.list(pvals_time[,2]), ])

head(pvals_time[sort.list(pvals_time[,3]), ])

gene <- "KCNN1"
gene <- "AP2S1"
expr_gene <- expr[gene,]
aframe <- data.frame(expr_gene, meta)
tmp <- aframe[which(aframe$patient %in% c("V", "C")),]
ggplot(tmp, aes(Days.of.Ibrutinib.Treatment, expr_gene, color = Ibrutinib.sensitivity)) +
  #facet_wrap(~ patient) +
  geom_smooth() +
  theme_bw() + ggtitle(gene)

# Run pseudobulk level analysis ####
asplit <- split(rownames(meta), meta$batch)
good <- names(which(unlist(lapply(asplit, length)) > 10))
asplit <- asplit[good]
pseudobulks <- do.call(cbind, lapply(asplit, function(x) rowMeans(expr[,x])))

meta_pseudo <- meta[match(colnames(pseudobulks), meta$batch),]

# Run model ####
pvals_time <- do.call(rbind, lapply(rownames(pseudobulks), function(gene){
  expr_gene <- pseudobulks[gene,]
  aframe <- data.frame(expr_gene, meta_pseudo)
  
  afit <- try(gam(expr_gene ~ Days.of.Ibrutinib.Treatment * Ibrutinib.sensitivity, data = aframe))
  
  if(class(afit) == "try-error") return(c(NA, NA, NA))
  summary(afit)[[4]][-1]
}))
rownames(pvals_time) <- rownames(pseudobulks)

# Visualize genes ####
head(pvals_time[sort.list(pvals_time[,1]), ])
gene <- "CCL3L1"
expr_gene <- pseudobulks[gene,]
aframe <- data.frame(expr_gene, meta_pseudo)
ggplot(aframe, aes(Days.of.Ibrutinib.Treatment, expr_gene, color = batch)) +
  #facet_wrap(~ Ibrutinib.sensitivity) +
  geom_point() +
  theme_bw() + ggtitle(gene)

head(pvals_time[sort.list(pvals_time[,2]), ])
gene <- "GPR82"
expr_gene <- pseudobulks[gene,]
aframe <- data.frame(expr_gene, meta_pseudo)
ggplot(aframe, aes(Ibrutinib.sensitivity, expr_gene, color = Ibrutinib.sensitivity)) +
  geom_boxplot() +
  theme_bw() + ggtitle(gene)

head(pvals_time[sort.list(pvals_time[,3]), ])
gene <- "NKG7"
expr_gene <- pseudobulks[gene,]
aframe <- data.frame(expr_gene, meta_pseudo)
ggplot(aframe, aes(Days.of.Ibrutinib.Treatment, expr_gene, color = patient)) +
  facet_wrap(~ Ibrutinib.sensitivity) +
  geom_point() +
  xlim(-10, 500) +
  theme_bw() + ggtitle(gene)

# Run PCA ####
pca <- prcomp(t(pseudobulks))
udata <- umap::umap(pca$x, n_neighbors = 5)
aframe <- data.frame(meta_pseudo, pca$x, udata$layout)
ggplot(aframe, aes(PC1, PC2, color = chemistry, shape = Ibrutinib.sensitivity)) +
  geom_point() +
  theme_bw()

ggplot(aframe, aes(PC3, PC4, color = patient, shape = Ibrutinib.sensitivity)) +
  geom_point() +
  theme_bw()

ggplot(aframe, aes(X1, X2, color = chemistry, shape = Ibrutinib.sensitivity)) +
  geom_point() +
  theme_bw()

# Supervised analysis using LDA ####
library(MASS)
aframe <- data.frame(outcome = meta_pseudo$Ibrutinib.sensitivity, t(pseudobulks))

lda_fit <- lda(outcome~., data = aframe[aframe$outcome %in% c("R", "S"),])

coords <- predict(lda_fit, aframe)
aframe <- data.frame(meta_pseudo, lda = coords$x, index = 1:nrow(meta_pseudo))
ggplot(aframe, aes(index, LD1, color = patient, shape = Ibrutinib.sensitivity)) +
  geom_point() +
  theme_bw()

# Try Svensson code ####
library(lmerTest)

# Load data ####
seu <- readRDS("/Users/lukas.simon/OneDrive/Miko/UTHealth/projects/Wang_collab/b_cell.rds")
meta <- seu@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1", "Normal_2"))] <- "no"
meta$batch <- as.factor(as.character(meta$batch))

f1 <- 'expr ~ offset(log(total_counts)) + (1 | patient)'
f2 <- 'expr ~ offset(log(total_counts)) + (1 | tumor/patient)'

# Function for fitting the two models to each gene and doing an LR test.
fit_lme <- function(gene = "CCND1") {
  raw_counts <- seu@assays$RNA@counts[gene,]
  tmp_df <- cbind(meta, data.frame(expr = raw_counts))
  
  mres1 <- glmer(f1, data = tmp_df, family = poisson())
  mres2 <- glmer(f2, data = tmp_df, family = poisson())
  anova_res <- anova(mres1, mres2)
  pval <- anova_res$`Pr(>Chisq)`[2]
  effect_size <- ranef(mres2)$tumor[,1][2] - ranef(mres2)$tumor[,1][1]
  
  c(pval = pval, effect_size = effect_size)
}

e_fit_lme <- function(gene) {
  afit <- try(fit_lme(gene))
  if(class(afit) == "try-error") return(c(pval = NA, effect_size = NA))
  afit
}

# Apply the function
genes <- names(which(rowSums(seu@assays$RNA@counts) > 2000))
genes <- sample(genes, 100)
fit_results <- do.call(rbind, lapply(genes, e_fit_lme))
rownames(fit_results) <- genes

# Visualize results
gene <- "CDC42"
expr_gene <- seu@assays$SCT@data[gene,]
aframe <- data.frame(expr_gene, meta)
ggplot(aframe, aes(tumor, expr_gene, color = patient)) +
  geom_boxplot() + theme_bw()

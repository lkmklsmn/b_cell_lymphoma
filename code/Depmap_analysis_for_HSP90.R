# Load R libs ####
library(data.table)
library(ggplot2)

# Load data ####
sample_info <- read.csv("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/sample_info.csv")
crispr <- fread("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/CRISPR_gene_effect.csv", sep = ",")
scores <- data.matrix(crispr[,-1])
rownames(scores) <- crispr[[1]]
genes <- colnames(scores)
genes <- unlist(lapply(genes, function(x) strsplit(x, " ",  fixed = T)[[1]][1]))
colnames(scores) <- genes

# Find cell lines present in both data tables ####
ok <- intersect(sample_info$DepMap_ID, rownames(scores))

# Calculate correlation between MYC and HSP90AA1/HSP90AB1 ####
asplit <- split(1:nrow(scores), sample_info[match(ok, sample_info$DepMap_ID), 'primary_disease'])
good_tissues <- names(which(unlist(lapply(asplit, length)) > 25))
tmp <- lapply(good_tissues, function(x){
  cor(scores[asplit[[x]], ], scores[asplit[[x]], "MYC"])
})

tmp <- do.call(rbind, lapply(tmp, function(x) rank(x[,1])[genes]))
rownames(tmp) <- good_tissues

aframe <- data.frame(primary_disease = good_tissues,
                     tmp)
aframe$type <- "NA"
aframe$type[aframe$primary_disease == "Lymphoma"] <- "MCL"
  
aframe$primary_disease <- factor(aframe$primary_disease,
                                 levels = aframe$primary_disease[order(aframe$HSP90AB1)])
p_hsp90ab1 <- ggplot(aframe, aes(HSP90AB1, primary_disease, fill = type)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c("red", "grey")) +
  theme_bw()

aframe$primary_disease <- factor(aframe$primary_disease,
                                 levels = aframe$primary_disease[order(aframe$HSP90AA1)])
p_hsp90aa1 <- ggplot(aframe, aes(HSP90AA1, primary_disease, fill = type)) +
  geom_bar(stat = 'identity') +
  scale_fill_manual(values = c("red", "grey")) +
  theme_bw()

gridExtra::grid.arrange(p_hsp90ab1, p_hsp90aa1, ncol = 2)

# ####
aframe <- data.frame(myc = scores[, "MYC"],
                     hsp90b1 = scores[, "HSP90B1"],
                     hsp90aa1 = scores[, "HSP90AA1"],
                     hsp90ab1 = scores[, "HSP90AB1"],
                     sample_info[match(ok, sample_info$DepMap_ID),])

p_hsp90ab1 <- ggplot(aframe[aframe$primary_disease %in% good_tissues,], aes(myc, hsp90ab1)) +
  facet_wrap(~ primary_disease) +
  geom_point() +
  xlab("MYC dependency") + ylab("HSP90AB1 dependency") +
  theme_bw()

p_hsp90aa1 <- ggplot(aframe[aframe$primary_disease %in% good_tissues,], aes(myc, hsp90aa1)) +
  facet_wrap(~ primary_disease) +
  geom_point() +
  xlab("MYC dependency") + ylab("HSP90AA1 dependency") +
  theme_bw()

gridExtra::grid.arrange(p_hsp90ab1, p_hsp90aa1, ncol = 2)


ggplot(aframe[aframe$primary_disease == "Lymphoma",], aes(myc, hsp90ab1)) +
  facet_wrap(~ primary_disease) +
  geom_smooth(method=lm) + geom_point() +
  xlab("MYC dependency") + ylab("HSP90AB1 dependency") +
  theme_bw()

ggplot(aframe[aframe$primary_disease == "Lymphoma",], aes(myc, hsp90aa1)) +
  facet_wrap(~ primary_disease) +
  geom_smooth(method=lm) + geom_point() +
  xlab("MYC dependency") + ylab("HSP90AB1 dependency") +
  theme_bw()

# ####
cell_lines <- sample_info$DepMap_ID[which(sample_info$primary_disease == "Lymphoma")]
cell_lines <- intersect(rownames(scores), cell_lines)

par(mfrow = c(1, 2))

correl <- cor(scores, scores[, "MYC"])
correl <- rev(sort(correl[,1]))[-1]
plot(density(correl, na.rm = T),
     main = "Across all cell lines",
     xlab = "Pearson correlation with MYC dependency")
abline(v = correl[c("HSP90AA1")], col = "red")
abline(v = correl[c("HSP90AB1")], col = "blue")

correl <- cor(scores[cell_lines, ], scores[cell_lines, "MYC"])
correl <- rev(sort(correl[,1]))[-1]
plot(density(correl, na.rm = T),
     main = "Across Lymphoma cell lines only",
     xlab = "Pearson correlation with MYC dependency")
abline(v = correl[c("HSP90AA1")], col = "red")
abline(v = correl[c("HSP90AB1")], col = "blue")
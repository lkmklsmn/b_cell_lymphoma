# Load packages ####
library(FNN)
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(readxl)
library(TCGAbiolinks)
library(tidyverse)

# Load global environment ####
load("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/global.RData")

# Correlation Matrices for pan and lym correlations with MYC ####
aframeCOR <- cor(kronos, kronos[, "MYC"])

ok <- intersect(rownames(sample_info[sample_info$lineage == "lymphocyte",]),
                rownames(kronos))
aframeLymCOR <- cor(kronos[ok, ], kronos[ok, "MYC"])

# Data frame with pan and lym correlations with MYC as columns ####
dfw <- as.data.frame(cbind(aframeCOR, aframeLymCOR))
names(dfw) <- c('pan_cancer', 'lymphocyte')

# Remove MYC row ####
dfw <- dfw[-which(row.names(dfw) %in% c('MYC')), ]
dfw <- rownames_to_column(dfw, var = "gene")

# Mininmil - Pan-cancer MYC Correlation Plot ####
dfw$gene <- factor(dfw$gene, levels = dfw$gene[order(-dfw$pan_cancer)])
dfw <- dfw[order(dfw$gene),]
p_pan <- ggplot(data = dfw, mapping = aes(x=gene, y=pan_cancer), color = NA) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_fill_manual(values = "lightgrey") +
  geom_hline(yintercept = 0, colour="black", linetype = 'dashed') +
  geom_vline(xintercept = which(dfw$gene == 'HSP90AB1'),
             colour = "red", size = 2) + 
  labs(title = 'Pan-cancer (n=1048)',
       x = 'Sorted genes',
       y = "Correlation with MYC") +
  theme_minimal()

# Mininmil - Lymphoma MYC Correlation Plot ####
dfw$gene <- factor(dfw$gene, levels = dfw$gene[order(-dfw$lymphocyte)])
dfw <- dfw[order(dfw$gene), ]
p_lym <- ggplot(data = dfw, mapping = aes(x=gene, y=lymphocyte), color = NA) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_fill_manual(values = "lightgrey") +
  geom_hline(yintercept = 0, colour="black", linetype = 'dashed') +
  geom_vline(xintercept = which(dfw$gene == 'HSP90AB1'),
             colour = "red", size = 2) + 
  labs(title = 'Lymphoma (n=35)',
       x = 'Sorted genes',
       y = "Correlation with MYC") +
  theme_minimal()

# Tertile plot ####
plot_tertiles_tissue <- function(gene = "HSP90AB1", tissue = "Lymphoma", refgene = "MYC"){
  aframe <- data.frame(gene = kronos[, gene],
                       refgene = kronos[, refgene],
                       sample_info[rownames(kronos),])
  
  good_tissues <- names(which(table(aframe$lineage) > 20))
  aframe <- aframe[aframe$lineage %in% good_tissues, ]
  
  if(tissue != "pan-cancer" & tissue != "all"){
    aframe <- aframe[aframe$primary_disease == tissue,]
  }
  if(tissue == "pan-cancer"){
    aframe$lineage <- "Pan-cancer"
  }
  
  aframe$split <- cut(aframe[, "refgene"],
                      breaks=c(quantile(aframe[, "refgene"],
                                        probs = seq(0, 1, length = 4))),
                      labels = c("low", "int", "high"), include.lowest = TRUE)
  aframe$split <- gsub("low", "dependent", aframe$split)
  aframe$split <- gsub("high", "independent", aframe$split)
  asplit <- split(aframe$gene, aframe$split)
  
  ggplot(aframe[aframe$split != "int", ],
         aes(split, gene, color = split)) +
    facet_wrap(~primary_disease) +
    geom_violin() +
    geom_jitter(width = .1) + 
    stat_compare_means(method =  "t.test") + # adds t-test p-values to plots
    scale_color_manual(values = c("red", "black"), name = "Status") +
    labs(y = gene, x = paste(refgene, "dependency")) +
    theme_classic()
}
p_box <- plot_tertiles_tissue(gene = "HSP90AB1", tissue = "Lymphoma", refgene = "MYC")

# Compile into single plot ####
grid.arrange(p_pan, p_lym, p_box, ncol = 3)

# Highlight specific cell lines in barplots ####
aframe <- data.frame(kronos[, c("MYC", "HSP90AB1")],
                     sample_info[rownames(kronos),])
aframe <- aframe[aframe$primary_disease == "Lymphoma", ]

aframe$split <- cut(aframe[, "MYC"],
                    breaks=c(quantile(aframe[, "MYC"],
                                      probs = seq(0, 1, length = 4))),
                    labels = c("low", "int", "high"), include.lowest = TRUE)
aframe$split <- gsub("low", "dependent", aframe$split)
aframe$split <- gsub("high", "independent", aframe$split)
aframe$split <- factor(aframe$split, levels = c("dependent", "int", "independent"))
aframe$cnv <- cnv[match(rownames(aframe), rownames(cnv)), "MYC"]

tmp <- aframe[, c("cnv", "MYC", "HSP90AB1", "split", "stripped_cell_line_name")]
tmp <- tmp[order(tmp$split, tmp$HSP90AB1), ]
head(tmp)
tail(tmp)

cell_lines <- c("C8166", "HDMYZ", "OCILY19", "A3KAW", "SMZ1", "RAJI")

aframe <- aframe[aframe$stripped_cell_line_name %in% cell_lines, ]
aframe$stripped_cell_line_name <- factor(aframe$stripped_cell_line_name,
                                         levels = cell_lines)

ggplot(aframe) +
  labs(y = "HSP90AB1 dependency",
       x = "Cell lines") +
  geom_bar(aes(x = stripped_cell_line_name, y = HSP90AB1, fill = split),
           stat = "identity") +
  scale_fill_manual(values = c("red", "black"), name = "MYC") +
  theme_classic() +
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5))
ggsave(
  filename = "/Users/lukas/OneDrive/Documents/GitHub/b_cell_lymphoma/figs/HSP90AB1_cell_lines.pdf",
  width = 4, height = 4)


# Download TCGA DLBC expression ####
query <- GDCquery(
  project = "TCGA-DLBC",
  data.category = "Gene expression",
  data.type = "Gene expression quantification",
  file.type  = "normalized_results",
  experimental.strategy = "RNA-Seq",
  legacy = TRUE
)
GDCdownload(
  query = query, 
  method = "api", 
  files.per.chunk = 10
)
data <- GDCprepare(query = query)

# Load cellinger data ####
celllinger <- read_excel("/Users/lukas/OneDrive/Documents/GitHub/b_cell_lymphoma/data/41467_2020_20294_MOESM4_ESM.xlsx", sheet = 1)
batch <- celllinger[["sampleID"]]
batch <- unlist(lapply(batch, function(x) strsplit(x, "-", fixed = T)[[1]][1]))
batch <- unlist(lapply(batch, function(x) strsplit(x, "_", fixed = T)[[1]][1]))
celllinger$batch <- batch

asplit <- split(1:nrow(coord), celllinger$lineage)
median_coord <- do.call(rbind,
                        lapply(asplit, function(x) c(UMAP_1 = median(coord[x, 1]),
                                                     UMAP_2 = median(coord[x, 2]))))
median_coord <- data.frame(label = names(asplit), median_coord)

coord <- celllinger[, c("UMAP_1", "UMAP_2")]
coord <- data.matrix(coord)

p_global <- ggplot(celllinger, 
       aes(UMAP_1, UMAP_2, fill = lineage,
           color = type)) +
  geom_point(pch = 21, alpha = 0.7)  +
  theme_classic() + 
  theme(legend.position = 'bottom', 
        legend.margin = margin(1,1,1,1)) +
  guides(fill = "none") +
  scale_color_manual(values = c("black", "white")) +
  labs(title = "Cellinger: global embedding",
       subtitle = "Maps cell lines to tumor samples",
       x = "UMAP 1", y = "UMAP 2") +
  ggrepel::geom_label_repel(data = median_coord, aes(label = label),
             fill = NA,
             color = 'black',
             size = 4)

kdata <- knn.index(coord, k = 20)
tcga_samples <- which(batch == "TCGA")
ok <- which(celllinger$batch == "ACH" & celllinger$disease == "Lymphoma")

neighbors <- apply(kdata[ok, ], 1 , function(x) intersect(x, tcga_samples))
neighbors <- unlist(neighbors)

tmp <- celllinger[union(ok, neighbors),]
tmp$group <- "Depmap"
tmp$group[tmp$batch == "TCGA"] <- tmp$subtype[tmp$batch == "TCGA"]
tmp$group[-which(tmp$group %in%
                   c("Depmap", "diffuse large B-cell lymphoma"))] <- "TCGA-Other"
tmp$group <- gsub("diffuse large B-cell lymphoma", "TCGA- DLBC",
                  tmp$group)

p_local <- ggplot(tmp,
       aes(UMAP_1, UMAP_2,
           color = group,
           shape = type)) +
  geom_point(alpha = 0.7, size = 3) +
  labs(title = "Celllinger mapping of Lymphoma cells",
       subtitle = "") +
  xlim(-2, 2) +
  ylim(-1.5, 2) +
  scale_color_manual(values = c("grey", "red", "orange"),
                    name = c("Depmap", "TCGA-DLBC", "TCGA-Other")) +
  theme_classic() +
  theme(legend.position = 'bottom', 
        legend.margin = margin(1,1,1,1))

# Plot MYC HSP90AB1 correlation in  TCGA ####
expr <- assays(data)[["normalized_count"]]
aframe <- data.frame(t(expr[c("HSP90AB1", "MYC"),]))
p_cor <- ggplot(aframe, aes(log(MYC), log(HSP90AB1))) +
  geom_point() + geom_smooth(method = lm) +
  labs(title = "TCGA-DLBC",
       subtitle = "RNA-seq",
       x = "MYC levels [log]",
       y = "HSP90AB1 levels [log]") +
  stat_cor() +
  theme_bw()  

# Compile into single plot ####
grid.arrange(p_global, p_local, p_cor, ncol = 3)

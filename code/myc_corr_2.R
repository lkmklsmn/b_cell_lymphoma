# Load global environment ####
load("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/global.RData")

# Load packages ####
library(gridExtra)
library(ggplot2)
library(ggpubr)
library(tidyverse)

################################################################################
# aframe Data Frame
aframe <- data.frame(kronos,
                     sample_info[rownames(kronos), ])

# Separate data frame filtered for lymphocyte tissues
aframeLym <- aframe %>% filter(lineage=='lymphocyte')

# Correlation Matrices for pan and lym correlations with MYC
aframeCOR <- aframe[,1:dim(kronos)[2]] %>% cor()
aframeCOR <- aframeCOR[, 'MYC']

aframeLymCOR <- aframeLym[,1:dim(kronos)[2]] %>% cor()
aframeLymCOR <- aframeLymCOR[, 'MYC']

################################################################################
# Data frame with pan and lym correlations with MYC as columns
dfw <- as.data.frame(cbind(aframeCOR,aframeLymCOR))
names(dfw) <- c('pan_cancer','lymphocyte')

# Remove MYC row
dfw <- dfw[!(row.names(dfw) %in% c('MYC')), ]
dfw <- rownames_to_column(dfw, var = "gene")

################################################################################
# Mininmil - Pan-cancer MYC Correlation Plot

dfw$gene <- factor(dfw$gene, levels = dfw$gene[order(-dfw$pan_cancer)])
dfw <- dfw[order(dfw$gene),]
p_pan <- ggplot(data = dfw, mapping = aes(x=gene, y=pan_cancer), color = NA) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_fill_manual(values = "lightgrey") +
  geom_hline(yintercept = 0, colour="black", linetype = 'dashed') +
  geom_vline(xintercept=which(dfw$gene == 'HSP90AB1'), colour="red", size = 2) + 
  labs(title='Pan-cancer (n=1048)', x='Sorted genes',y="Correlation with MYC") +
  theme_minimal()

################################################################################
# Mininmil - Lymphoma MYC Correlation Plot

dfw$gene <- factor(dfw$gene, levels = dfw$gene[order(-dfw$lymphocyte)])
dfw <- dfw[order(dfw$gene),]
p_lym <- ggplot(data = dfw, mapping = aes(x=gene, y=lymphocyte), color = NA) +
  geom_bar(stat = "identity") +
  scale_x_discrete(labels = NULL, breaks = NULL) +
  scale_fill_manual(values = "lightgrey") +
  geom_hline(yintercept = 0, colour="black", linetype = 'dashed') +
  geom_vline(xintercept = which(dfw$gene == 'HSP90AB1'), colour = "red", size = 2) + 
  labs(title='Lymphoma (n=35)', x='Sorted genes',y="Correlation with MYC") +
  theme_minimal()

################################################################################
#  Tertile plot 
plot_tertiles_tissue <- function(gene = "HSP90AB1", tissue = "Lymphoma", refgene = "MYC"){
  aframe <- data.frame(gene = kronos[, gene],
                       refgene = kronos[, refgene],
                       sample_info[rownames(kronos),])
  
  good_tissues <- names(which(table(aframe$lineage) > 20))
  aframe <- aframe[aframe$lineage %in% good_tissues,]
  
  if(tissue != "pan-cancer" & tissue != "all"){
    aframe <- aframe[aframe$primary_disease == tissue,]
  }
  if(tissue == "pan-cancer"){
    aframe$lineage <- "Pan-cancer"
  }
  
  aframe$split <- cut(aframe[, "refgene"], breaks=c(quantile(aframe[, "refgene"],
                                                             probs = seq(0, 1, length = 4))),
                      labels = c("low", "int", "high"), include.lowest = TRUE)
  aframe$split <- gsub("low", "dependent", aframe$split)
  aframe$split <- gsub("high", "independent", aframe$split)
  asplit <- split(aframe$gene, aframe$split)
  
  ggplot(aframe[aframe$split != "int", ], aes(split, gene, color = split)) +
    facet_wrap(~primary_disease) +
    geom_violin() +
    geom_jitter(width = .1) + 
    stat_compare_means(method =  "t.test") + # adds t-test p-values to plots
    scale_color_manual(values = c("red", "black"), name = "Status") +
    labs(y = gene, x = paste(refgene, "dependency")) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) # adds space above plots
}
p_box <- plot_tertiles_tissue(gene = "HSP90AB1", tissue = "Lymphoma", refgene = "MYC")

################################################################################
grid.arrange(p_pan, p_lym, p_box, ncol = 3)

################################################################################
mut <- data.table::fread("/Users/lukas/OneDrive/Documents/GitHub/depmap_app/data/CCLE_mutations.csv", sep = ",")
snps <- unique(mut$Genome_Change[grep("HSP90AB1", mut$Hugo_Symbol)])
cells <- mut$DepMap_ID[which(mut$Genome_Change %in% snps)]
aframe <- data.frame(kronos[, c("MYC", "HSP90AB1")],
                     sample_info[rownames(kronos),])
aframe$mutant <- 'no'
aframe$mutant[aframe$DepMap_ID %in% cells] <- "yes"
ggplot(aframe,
       aes(mutant, MYC, color = mutant)) +
  geom_violin() + geom_point() +
  theme_bw()

aframe <- data.frame(rnaseq[, c("MYC", "HSP90AB1", "HSP90AA1", "HSP90B1")],
                     sample_info[rownames(rnaseq),])
ggplot(aframe[aframe$lineage == 'lymphocyte',],
       aes(MYC, HSP90B1)) +
  geom_smooth(method = lm) + geom_point() + 
  stat_cor() +
  theme_bw()

aframe <- data.frame(cnv[, c("HSP90AB1", "HSP90AA1", "HSP90B1")],
                     MYC = kronos[match(rownames(cnv), rownames(kronos)), "MYC"],
                     sample_info[rownames(cnv),])
ggplot(aframe[aframe$lineage == "lymphocyte",],
       aes(HSP90AB1 > 1.2, MYC)) +
  geom_violin() + geom_point() + 
  stat_compare_means(method = "t.test") +
  theme_bw()

aframe <- data.frame(rnaseq[, c("MYC", "HSP90AB1", "HSP90AA1", "HSP90B1")],
                     kronos[match(rownames(rnaseq), rownames(kronos)), c("MYC", "HSP90AB1", "HSP90AA1", "HSP90B1")],
                     sample_info[rownames(rnaseq),])
ggplot(aframe[aframe$lineage == "lymphocyte",],
       aes(HSP90AB1.1, MYC.1, color = HSP90AB1 + MYC)) +
  geom_point() + 
  stat_cor() +
  theme_bw()

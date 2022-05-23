setwd("/Users/yanfangfang/Downloads/MW/")
# Load R libs ####
library(data.table)
library(ggplot2)
library(patchwork)
library(ggpubr)

# Load data ####
sample_info <- read.csv("data/depmap/sample_info.csv")
cn <- fread("data/depmap/CRISPR_gene_effect.csv", sep = ",")
scores <- data.matrix(cn[,-1])
rownames(scores) <- cn$DepMap_ID
genes <- colnames(scores)
genes <- unlist(lapply(genes, function(x) strsplit(x, " ",  fixed = T)[[1]][1]))
colnames(scores) <- genes

# Find cell lines present in both data tables ####
ok <- intersect(sample_info$DepMap_ID, rownames(scores))
scores <- scores[ok,]
sample_info <- sample_info[match(ok, sample_info$DepMap_ID),]

# tmp <- sample_info[sample_info$primary_disease == "Lymphoma",]
# tmp <- tmp[tmp$lineage_sub_subtype=='b_cell_mantle_cell',]

# Calculate correlation between MYC and HSP90AA1/HSP90AB1 ####
gene1 <- "MYC"
gene2 <- "HSP90AB1"

# Across all primary disease ####
asplit <- split(1:nrow(scores), sample_info[, 'primary_disease'])
good_tissues <- names(which(unlist(lapply(asplit, length)) > 25))
aframe <- data.frame('gene1' = scores[, gene1],'gene2' = scores[, gene2],sample_info)
ggplot(aframe[aframe$primary_disease %in% good_tissues,], aes(gene1, gene2)) +
        geom_point() +geom_smooth(method='lm')+
        xlab(paste0(gene1," dependency")) + ylab(paste0(gene2," dependency"))+theme_bw()+
        stat_cor(method = "pearson")+ggtitle("Across all cell lines")
p1 <- ggplot(aframe[aframe$primary_disease %in% good_tissues,], aes(gene1, gene2)) +
        facet_wrap(~ primary_disease,scales = 'free') +
        geom_point() +geom_smooth(method='lm')+
        xlab(paste0(gene1," dependency")) + ylab(paste0(gene2," dependency")) +
        theme_bw()+ggtitle("Across all cell lines")
library(ggrepel)
tmp <- aframe[aframe$primary_disease == "Lymphoma" & aframe$lineage_sub_subtype!="",]
tmp$MCL <- "other"
tmp$MCL[tmp$lineage_sub_subtype=="b_cell_mantle_cell"] <- "MCL"
p2 <- ggplot(tmp, aes(gene1, gene2)) +
        geom_smooth() + geom_point(aes(color=MCL,shape=lineage_sub_subtype)) +
        geom_label_repel(aes(label=lineage_sub_subtype),size=2,max.overlaps = 3)+
        xlab(paste0(gene1," dependency")) + ylab(paste0(gene2," dependency")) +
        theme_bw()+
        stat_cor(method = "pearson")+
        ggtitle("Lymphoma cell lines only")

# bar plot ####
tmp <- lapply(good_tissues, function(x){
        cor(scores[asplit[[x]], ], scores[asplit[[x]], gene1])
})
tmp <- do.call(rbind, lapply(tmp, function(x) rank(x[,1])[genes]))
rownames(tmp) <- good_tissues
aframe <- data.frame('primary_disease' = good_tissues,tmp)
aframe$type <- "others"
aframe$type[aframe$primary_disease == "Lymphoma"] <- "Lymphoma"
aframe$primary_disease <- factor(aframe$primary_disease,
                                 levels = aframe$primary_disease[order(aframe[,gene2])])
aframe$tmp <- aframe[,gene2]
p3 <- ggplot(aframe, aes(tmp, primary_disease, fill = type)) +
        geom_bar(stat = 'identity') + scale_fill_manual(values = c("red", "grey")) +
        theme_bw()+ylab("Primary disease")+
        xlab(paste0("correlation of ",gene2," with ",gene1))+
        ggtitle("Across all cell lines")

# density plot ####
correl <- cor(scores, scores[, gene1])
correl <- rev(sort(correl[,1]))[-1]
df <- data.frame("corr"=correl)
intercept <- df[gene2,]
p4 <- ggplot(df, aes(x=corr)) + geom_density()+
        geom_vline(xintercept=intercept,color="steelblue", linetype="dashed", size=1)+
        ggtitle("Across all cell lines")+
        xlab(paste0("Pearson correlation with ",gene1," dependency"))+
        ylab("Density")+theme_bw()

cell_lines <- sample_info$DepMap_ID[which(sample_info$primary_disease == "Lymphoma")]
cell_lines <- intersect(rownames(scores), cell_lines)
correl <- cor(scores[cell_lines, ], scores[cell_lines, gene1])
correl <- rev(sort(correl[,1]))[-1]
df <- data.frame("corr"=correl) 
intercept <- df[gene2,]
p5 <- ggplot(df, aes(x=corr)) + geom_density()+
        geom_vline(xintercept=intercept,color="steelblue", linetype="dashed", size=1)+
        ggtitle("Across lymphoma cell lines only")+
        xlab(paste0("Pearson correlation with ",gene1," dependency"))+
        ylab("Density")+theme_bw()
# p <- ggplot_build(p5)
# h_value <- max(p$data[[1]]$density)-0.2
# p5 <- p5+geom_text(aes(intercept,h_value,label=gene2))
p4|p5

# Within lymphoma ####
sub_sample_info <- sample_info[sample_info$primary_disease=="Lymphoma" &
                                       sample_info$Subtype!="",]
ok <- intersect(sub_sample_info$DepMap_ID, rownames(scores))
sub_scores <- scores[ok,]
asplit <- split(1:nrow(sub_scores), sub_sample_info[match(ok, sub_sample_info$DepMap_ID), 
                                                    'Subtype'])
good_tissues <- names(which(unlist(lapply(asplit, length)) > 0))
tmp <- lapply(good_tissues, function(x){
        cor(sub_scores[asplit[[x]], ,drop=F], sub_scores[asplit[[x]], "MYC",drop=F])
})
tmp <- do.call(rbind, lapply(tmp, function(x) rank(x[,1])[genes]))
rownames(tmp) <- good_tissues
aframe <- data.frame(subtype = good_tissues,tmp)
aframe$subtype <- factor(aframe$subtype,
                         levels = aframe$subtype[order(aframe$HSP90AB1)])
p6 <- ggplot(aframe, aes(HSP90AB1, subtype)) +
        geom_bar(stat = 'identity') + scale_fill_manual(values = c("red", "grey")) +
        theme_bw()+ylab("Lymphoma subtypes")+
        ggtitle("Across subtypes of lymphoma cell lines")+
        theme(legend.position = 'none')

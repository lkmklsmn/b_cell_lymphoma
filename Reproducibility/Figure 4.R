setwd("/Users/yanfangfang/Downloads/MW/")

# Load R libs ####
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(zoo)

# Load Seurat object ####
load("../data/BlueYellowColormaps_V1.RData")
b_seu <- readRDS("../data/integrated_b_cells.rds")
good <- names(table(b_seu$sample))[(table(b_seu$sample) > 20)]
b_seu <- b_seu[,b_seu$sample %in% good]
b_seu$clinical_outcome <- factor(b_seu$ibrutinib_sensitivity,levels = c("Normal","S","Slow_responder","R","Dual"))
levels(b_seu$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")

# Figure 4A, Outcome-specific DEGs ####
res <- read.delim("../outputs/mixed_model_DEG_results.txt",sep="\t",row.names = 1)
df <- na.omit(res)
df$dual_spe <- apply(df,1,function(x){x[5]-mean(c(x[6],x[7],x[8]))})
df$r_spe <- apply(df,1,function(x){x[6]-mean(c(x[5],x[7],x[8]))})
df$s_spe <- apply(df,1,function(x){x[7]-mean(c(x[5],x[6],x[8]))})
df$slow_spe <- apply(df,1,function(x){x[8]-mean(c(x[5],x[6],x[7]))})
no <- 100
# dual-specific
tmp <- df[df$coef_Dual>0.2,]
tmp <- tmp[order(tmp$dual_spe,decreasing = T),]
dual <-  rownames(tmp)[1:no]
# r-specific
tmp <- df[df$coef_R>0.2,]
tmp <- tmp[order(tmp$r_spe,decreasing = T),]
r <-  rownames(tmp)[1:no]
# slow-specific
tmp <- df[df$coef_slow>0.2,]
tmp <- tmp[order(tmp$slow_spe,decreasing = T),]
slow <-  rownames(tmp)[1:no]
# s-specific
tmp <- df[df$coef_S>0.2,]
tmp <- tmp[order(tmp$s_spe,decreasing = T),]
s <-  rownames(tmp)[1:no]
s <- s[!(s %in% slow)]
de_genes <- c(s,slow,r,dual)

# Heatmap ####
cohort_1 <- b_seu[, b_seu$chemistry=="3prime"]
cohort_1 <- ScaleData(cohort_1,features = de_genes)
asplit_cells <- split(rownames(cohort_1@meta.data), cohort_1$sample)
len <- unlist(lapply(asplit_cells, length))
good <- names(len)[len > 50]
asplit_cells <- asplit_cells[good]
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(cohort_1@assays$RNA@scale.data[de_genes, sample(unlist(x), 20)])
        s2 <- Matrix::rowMeans(cohort_1@assays$RNA@scale.data[de_genes, sample(unlist(x), 20)])
        s3 <- Matrix::rowMeans(cohort_1@assays$RNA@scale.data[de_genes, sample(unlist(x), 20)])
        cbind(s1, s2, s3)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 3)))
anno_col <- data.frame(sample,
                       "Cohort" = cohort_1@meta.data[match(sample,cohort_1$sample),
                                                     'chemistry'],
                       "clinical_outcome" = cohort_1@meta.data[match(sample,cohort_1$sample),
                                                               'clinical_outcome'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col <- arrange(anno_col,factor(clinical_outcome,
                                    levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")))
lst <- unique(anno_col$sample)
anno_col$sample <- factor(anno_col$sample,levels=lst)
anno_col$Cohort <- "BTKi cohort"
annotation_colors = list(
        clinical_outcome = c('Normal' = "darkgrey",
                             'BTKi-Fast'="steelblue",
                             'BTKi-Slow'="lightblue2",
                             'BTKi-R'="peachpuff4",
                             'Dual-R'="chocolate2"),
        Cohort = c('BTKi cohort'="#F8766D"),
        DEGs = c('BTKi-Fast-specific'="steelblue",
                 'BTKi-Slow-specific'="lightblue2",
                 'BTKi-R-specific'="peachpuff4",
                 'Dual-R-specific'="chocolate2")
)
annotation_row <- data.frame("DEGs" = c(rep("BTKi-Fast-specific",length(s)),
                                        rep("BTKi-Slow-specific",length(slow)),
                                        rep("BTKi-R-specific",length(r)),
                                        rep("Dual-R-specific",length(dual))),
                             row.names = de_genes)

tmp <- do.call(cbind, lapply(split(rownames(annotation_row), annotation_row$DEGs), function(x) colMeans(means[x, ])))
tmp <- data.frame(tmp)
colnames(tmp) <- gsub(".specific", "", colnames(tmp), fixed = T)
colnames(tmp) <- gsub(".", "-", fixed = T, colnames(tmp))
tmp$score <- NA
lapply(rownames(tmp), function(x){
  print(x)
  outcome <- as.character(anno_col[x, "clinical_outcome"])
  try(tmp[x, "score"] <<- tmp[x, outcome])
})
tmp <- tmp[rownames(anno_col),]
tmp <- do.call(rbind, lapply(split(rownames(tmp), anno_col$sample), function(x) colMeans(tmp[x, ])))

anno_col$score <- tmp[anno_col$sample, "score"]
anno_col <- anno_col[order(anno_col$clinical_outcome, -anno_col$score),]

rownames(annotation_row) <- de_genes
means <- t(apply(means, 1, function(x) (x - mean(x))/sd(x)))
pheatmap(means[,rownames(anno_col)],cluster_rows = F, cluster_cols = F,
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         annotation_col = anno_col,
         annotation_row = annotation_row,
         show_colnames = F, show_rownames = F,
         annotation_colors = annotation_colors,
         border_color = NA,
         #gaps_row = c(length(s),length(s)+length(slow),length(s)+length(slow)+length(r)),
         main = "Outcome-specific differentially expressed genes (DEGs)") 

# Figure 4B ####
expr <- cohort_1@assays$RNA@data
meta <- cohort_1@meta.data

genes <- c("MYLIP", "FAM177B", "DDX11")
aframe <- data.frame(meta, t(expr[genes, ]))

asplit <- split(rownames(aframe), aframe$sample)
means <- do.call(cbind, lapply(asplit, function(x) rowMeans(expr[genes, x] > 0)))
aframe <- data.frame(aframe[match(colnames(means), aframe$sample), colnames(meta)],
                     t(means))

aframe <- reshape2::melt(aframe, measure.vars = genes)

meta_sample <- meta[match(unique(aframe$sample), meta$sample), ]
lst <- meta_sample$sample[order(meta_sample$clinical_outcome)]
aframe$sample <- factor(aframe$sample, levels = lst)

p_outcome <- ggplot(aframe, aes(clinical_outcome, value, fill = clinical_outcome)) +
  facet_wrap(~ variable, ncol = 1, scales = "free") + geom_boxplot() + geom_point() +
  scale_fill_manual(values = c("darkgrey", "steelblue", "lightblue2", "peachpuff4", "chocolate2"))+
  labs(y = "Expression levels",
       x = "Sample") +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.title.x = element_blank())+
  theme(legend.position = "bottom")

# Figure 4C ####
res_dual_r <- read.delim("../outputs/mixed_model_DEG_Dual_vs_R.txt", sep = "\t", row.names = 1)
res_dual_r <- na.omit(res_dual_r)
res_dual_r$adj.p <- p.adjust(res_dual_r$pval_Dual,method = "BH")
tmp <- res_dual_r[res_dual_r$adj.p < 0.1,]
tmp <- tmp[order(tmp$adj.p,decreasing = F),]
dual_r <-  rownames(tmp)
# heatmap
tmp <- b_seu[, b_seu$clinical_outcome %in% c("BTKi-R","Dual-R")]
good <- names(table(tmp$sample))[table(tmp$sample) > 50]
tmp <- tmp[, tmp$sample %in% good]
DefaultAssay(tmp) <- "RNA"
tmp$sample[tmp$sample == "A3" & tmp$chemistry == "5prime"] <- "A3_cohort2"
tmp$sample[tmp$sample == "A3" & tmp$chemistry == "3prime"] <- "A3_cohort1"
tmp <- ScaleData(tmp, features = rownames(tmp))
asplit_cells <- split(rownames(tmp@meta.data), tmp$sample)
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(tmp@assays$RNA@data[dual_r, sample(unlist(x), 50,replace=T)])
        s2 <- Matrix::rowMeans(tmp@assays$RNA@data[dual_r, sample(unlist(x), 50,replace=T)])
        s3 <- Matrix::rowMeans(tmp@assays$RNA@data[dual_r, sample(unlist(x), 50,replace=T)])
        cbind(s1, s2, s3)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 3)))
anno_col <- data.frame(sample,
                       "clinical_outcome" = tmp@meta.data[match(sample, tmp$sample),'clinical_outcome'],
                       "cohort" = tmp@meta.data[match(sample, tmp$sample),'chemistry'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col$cohort <- gsub("5prime", "CART cohort", anno_col$cohort)
anno_col$cohort <- gsub("3prime", "BTKi cohort", anno_col$cohort)
anno_col$cohort <- factor(anno_col$cohort, levels = c("BTKi cohort", "CART cohort"))
anno_col$clinical_outcome <- as.character(anno_col$clinical_outcome)
annotation_colors = list(
        clinical_outcome = c('BTKi-R' = "peachpuff4",
                             'Dual-R' = "chocolate2"),
        cohort = c('BTKi cohort' = "#F8766D", 'CART cohort' = "#00BFC4")
)


scale <- function(matr) t(apply(matr, 1, function(x) (x - mean(x))/sd(x)))
tmp <- cbind(scale(means[, anno_col$cohort == "BTKi cohort"]),
             scale(means[, anno_col$cohort == "CART cohort"]))
tmp <- tmp[which(apply(tmp, 1, var) > 0), ]

anno_col$score <- apply(tmp[, rownames(anno_col)], 2, median)
lvs <- names(sort(unlist(lapply(split(anno_col$score, anno_col$sample), mean))))
anno_col$sample <- factor(anno_col$sample, levels = lvs)
anno_col <- anno_col[order(anno_col$cohort, anno_col$clinical_outcome, anno_col$sample), ]

pheatmap(tmp[, rownames(anno_col)],
         cluster_rows = F, cluster_cols = F,
         breaks = seq(-2, 2, length = length(yellow2blue) + 1),
         col = yellow2blue, 
         gaps_col = length(which(anno_col$cohort == "BTKi cohort")),
         fontsize_row = 8,
         annotation_col = anno_col,
         show_colnames = F, show_rownames = T,
         annotation_colors = annotation_colors,
         main = "Dual vs IBN-R across both cohorts") 

# Figure 4D ####
create_plot_sample_level <- function(x){
  tmp <- b_seu[, b_seu$clinical_outcome %in% c("BTKi-R","Dual-R")]
  meta <- tmp@meta.data
  expr <- as.numeric(t(tmp@assays$RNA@data[x, ]))
  expr_mean <- unlist(lapply(split(expr, tmp$sample), mean))
  
  df <- data.frame(expr = expr_mean,
                   meta[match(names(expr_mean), meta$sample), ])
  df$Cohort <- gsub("3prime", "BTKi cohort", df$chemistry)
  df$Cohort <- gsub("5prime", "CART cohort", df$Cohort)
  df$clinical_outcome <- factor(df$clinical_outcome, levels = c("BTKi-R","Dual-R"))
  
  ggplot(df, aes(clinical_outcome, expr, fill = clinical_outcome)) +
    geom_boxplot() + geom_point() +
    scale_fill_manual(values = c("peachpuff4","chocolate2")) +
    ylab("Relative expression levels") +
    theme_bw()+xlab("Clinical outcome") +
    theme(axis.text.x = element_text(size=6, angle=45, h=1),
          axis.title.x = element_blank()) +
    ggtitle(x) +
    facet_wrap(~ Cohort, scales = "free")
}
p_dual <- grid.arrange(create_plot_sample_level('CDK9'),
                       create_plot_sample_level('POLR2C'))
grid.arrange(p_outcome, p_dual, ncol = 2)

# Figure 4E ####
df <- read.delim("mixed_model_DEG_GSEA/Project_Dual_R_both_cohorts/enrichment_results_Dual_R_both_cohorts.txt")
df <- df[df$FDR<0.05,c(1,4,6)]
colnames(df) <- c('geneSet','EnrichmentScore','FDR')
df$contrast <- "Dual-R vs BTKi-R"
df2 <- read.delim("mixed_model_DEG_GSEA/Project_R_vs_Slow_S/enrichment_results_R_vs_Slow_S.txt")
df2 <- df2[df2$FDR<0.05,c(1,4,6)]
colnames(df2) <- c('geneSet','EnrichmentScore','FDR')
df2$contrast <- "BTKi-R vs BTKi-S"
tmp <- data.frame(rbind(df,df2))
tmp$geneSet <- gsub("HALLMARK_","",tmp$geneSet)
library(tidytext)
library(dplyr)
tmp %>% mutate(geneSet = reorder_within(geneSet,EnrichmentScore, contrast)) %>%
        ggplot(aes(y=geneSet,x=EnrichmentScore))+geom_bar(stat="identity",fill='steelblue')+
        scale_y_reordered() +
        facet_wrap(~contrast,scales="free",ncol=1)+theme_bw()+ggtitle("Enriched pathways")+
        theme(axis.title.y =element_blank())

# Figure 4F ####
tmp <- read.delim("/Users/yanfangfang/h.all.v7.4.symbols.gmt",header = F,fill = NA)
rownames(tmp) <- tmp$V1
tmp <- tmp[,-c(1,2)]
x <- lapply(rownames(tmp),function(x){
        t <- unname(unlist(tmp[x,]))
        t[t!=""]
})
names(x) <- rownames(tmp)
oxphos <- x[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']]
myc_v1 <- x[['HALLMARK_MYC_TARGETS_V1']]
myc_v2 <- x[['HALLMARK_MYC_TARGETS_V2']]

seu <- b_seu
DefaultAssay(seu) <- "RNA"
seu$sample[seu$sample=="A3" & seu$chemistry=="3prime"] <- "A3_cohort1"
seu$sample[seu$sample=="A3" & seu$chemistry=="5prime"] <- "A3_cohort2"

myc_v1 <- myc_v1[myc_v1 %in% rownames(seu)]
myc_v2 <- myc_v2[myc_v2 %in% rownames(seu)]
oxphos <- oxphos[oxphos %in% rownames(seu)]

expr <- seu@assays$RNA@data
m <- mean(as.matrix(expr[oxphos,seu$clinical_outcome =="Normal"]))
OXPHOS <- colMeans(expr[oxphos,])/m
m <- mean(as.matrix(expr[myc_v1,seu$clinical_outcome=="Normal"]))
MYC_TARGETS_V1 <- colMeans(expr[myc_v1,])/m
m <- mean(as.matrix(expr[myc_v2,seu$clinical_outcome=="Normal"]))
MYC_TARGETS_V2 <- colMeans(expr[myc_v2,])/m
df1 <- data.frame(pathway="OXPHOS",expression=OXPHOS,
                  clinical_outcome=seu$clinical_outcome)
df2 <- data.frame(pathway="MYC_TARGETS_V1",expression=MYC_TARGETS_V1,
                  clinical_outcome=seu$clinical_outcome)
df3 <- data.frame(pathway="MYC_TARGETS_V2",expression=MYC_TARGETS_V2,
                  clinical_outcome=seu$clinical_outcome)
df <- rbind(df1,df2,df3)
ggplot(df, aes(x=factor(clinical_outcome,levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")), 
               y=expression,fill=pathway)) + 
        geom_boxplot(outlier.shape=NA)+ylim(0,6)+
        xlab("Clinical outcome")+ylab("Average pathway score")+
        ggtitle("Average expression of MYC and OXPHOS pathways")+
        theme_bw()+geom_hline(yintercept = 1,linetype="dashed")



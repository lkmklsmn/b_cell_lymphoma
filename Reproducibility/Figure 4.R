setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
library(pheatmap)
library(zoo)
load("/Users/yanfangfang/BlueYellowColormaps_V1.RData")
b_seu <- readRDS("data/integrated_b_cells.rds")
good <- names(table(b_seu$sample))[(table(b_seu$sample)>20)]
b_seu <- b_seu[,b_seu$sample %in% good]
b_seu$clinical_outcome <- factor(b_seu$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(b_seu$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")


# Figure 4A, Outcome-specific DEGs ####
res <- read.delim("mixed_model_DEG_results.txt",sep="\t",row.names = 1)
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
tmp <- b_seu[,b_seu$chemistry=="3prime"]
tmp <- ScaleData(tmp,features = de_genes)
asplit_cells <- split(rownames(tmp@meta.data), tmp$sample)
len <- unlist(lapply(asplit_cells,length))
good <- names(len)[len>50]
asplit_cells <- asplit_cells[good]
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[de_genes, sample(unlist(x), 20)])
        s2 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[de_genes, sample(unlist(x), 20)])
        s3 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[de_genes, sample(unlist(x), 20)])
        cbind(s1, s2, s3)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 3)))
anno_col <- data.frame(sample,"Cohort"=tmp@meta.data[match(sample,tmp$sample),'chemistry'],
                       "clinical_outcome"=tmp@meta.data[match(sample,tmp$sample),'clinical_outcome'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col <- arrange(anno_col,factor(clinical_outcome,
                                    levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")))
lst <- unique(anno_col$sample)
anno_col$sample <- factor(anno_col$sample,levels=lst)
anno_col$Cohort <- "BTKi cohort"
annotation_colors = list(
        clinical_outcome = c(Normal="darkgrey",'BTKi-Fast'="steelblue",
                              'BTKi-Slow'="lightblue2",'BTKi-R'="peachpuff4",
                              'Dual-R'="chocolate2"),
        Cohort=c('BTKi cohort'="#F8766D"),
        DEGs=c('BTKi-Fast-specific'="steelblue",'BTKi-Slow-specific'="lightblue2",
               'BTKi-R-specific'="peachpuff4",'Dual-R-specific'="chocolate2")
)
annotation_row <- data.frame("DEGs"=c(rep("BTKi-Fast-specific",length(s)),rep("BTKi-Slow-specific",length(slow)),
                                      rep("BTKi-R-specific",length(r)),rep("Dual-R-specific",length(dual))))
rownames(annotation_row) <- de_genes
pheatmap(means[,rownames(anno_col)],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         annotation_col = anno_col,show_colnames = F,annotation_colors=annotation_colors,
         show_rownames=F,
         annotation_row = annotation_row,
         #gaps_row = c(length(s),length(s)+length(slow),length(s)+length(slow)+length(r)),
         main="Outcome-specific differentially expressed genes (DEGs)") 

# Figure 4B ####
tmp <- ScaleData(tmp,features=rownames(tmp))
expr <- tmp@assays$RNA@scale.data
meta <- tmp@meta.data
create_plot <- function(gene = "CD52"){
        aframe <- data.frame(meta, expr = expr[gene,])
        aframe <- arrange(aframe,factor(clinical_outcome,
                                        levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")))
        lst <- unique(aframe$sample)
        ggplot(aframe, aes(factor(sample,levels=lst), expr, group = sample, fill = clinical_outcome)) +
                geom_boxplot(outlier.size = 0.2) +
                scale_fill_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
                ylab("Expression levels") +
                theme_bw()+xlab("sample") +
                theme(axis.text.x = element_text(size=8, angle=45),axis.title.x =element_blank(),
                      axis.title.y=element_blank())+
                theme(legend.position = "none")
}
create_plot(dual[4])+ggtitle(paste0("Dual-R-specific, ",dual[4])) ->p1
create_plot(r[6])+ggtitle(paste0("BTKi-R-specific, ",r[6])) ->p2
create_plot(slow[2])+ggtitle(paste0("BTKi-Slow-specific, ",slow[2]))->p3
create_plot(s[1])+ggtitle(paste0("BTKi-S-specific, ",s[1])) ->p4
p4/p3/p2/p1

# Figure 4C ####
res_dual_r <- read.delim("mixed_model_DEG_Dual_vs_R.txt",sep="\t",row.names = 1)
res_dual_r <- na.omit(res_dual_r)
res_dual_r$adj.p <- p.adjust(res_dual_r$pval_Dual,method="BH")
tmp <- res_dual_r[res_dual_r$adj.p<0.1,]
tmp <- tmp[order(tmp$adj.p,decreasing = F),]
dual_r <-  rownames(tmp)
# heatmap
tmp <- b_seu[,b_seu$clinical_outcome %in% c("BTKi-R","Dual-R")]
good <- names(table(tmp$sample))[table(tmp$sample)>50]
tmp <- tmp[,tmp$sample %in% good]
DefaultAssay(tmp) <- "RNA"
tmp$sample[tmp$sample=="A3" & tmp$chemistry=="5prime"] <- "A3_cohort2"
tmp$sample[tmp$sample=="A3" & tmp$chemistry=="3prime"] <- "A3_cohort1"
tmp <- ScaleData(tmp,features = rownames(tmp))
asplit_cells <- split(rownames(tmp@meta.data), tmp$sample)
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        s2 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        s3 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        s4 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        s5 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        cbind(s1, s2, s3,s4,s5)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 5)))
anno_col <- data.frame(sample,"clinical_outcome"=tmp@meta.data[match(sample,tmp$sample),'clinical_outcome'],
                       "cohort"=tmp@meta.data[match(sample,tmp$sample),'chemistry'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col$cohort <- as.factor(anno_col$cohort)
levels(anno_col$cohort) <- c("cohort1","cohort2")
anno_col$clinical_outcome <- as.character(anno_col$clinical_outcome)
anno_col <- arrange(anno_col,factor(clinical_outcome,levels=c("BTKi-R","Dual-R")))
annotation_colors = list(
        clinical_outcome = c('BTKi-R'="peachpuff4",'Dual-R'="chocolate2"),
        cohort=c(cohort1="#F8766D",cohort2="#00BFC4")
)
anno_col <- arrange(anno_col,factor(anno_col$cohort,levels=c("cohort1","cohort2")))
pheatmap(means[,rownames(anno_col)],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,gaps_col = 25,
         fontsize_row = 8,
         annotation_col = anno_col,show_colnames = F,show_rownames = T,annotation_colors=annotation_colors,
         main="Dual vs IBN-R across both cohorts") 

# Figure 4D ####
create_plot_cell_level <- function(gene = "CDK9"){
        tmp <- ScaleData(tmp,features = x)
        meta <- tmp@meta.data
        aframe <- data.frame(meta, expr = t(tmp@assays$RNA@scale.data))
        aframe <- arrange(aframe,factor(clinical_outcome,
                                        levels=c("BTKi-R","Dual-R")))
        lst <- unique(aframe$sample)
        aframe$chemistry <- as.factor(aframe$chemistry)
        levels(aframe$chemistry) <- c("BTKi cohort","CAR-T cohort")
        ggplot(aframe, aes(factor(sample,levels=lst), expr, group = sample, fill = clinical_outcome)) +
                geom_boxplot(outlier.size = 0.2) +
                scale_fill_manual(values=c("peachpuff4","chocolate2"))+
                ylab("Expression levels") +ggtitle(paste0(x," expression"))+
                theme_bw()+xlab("sample") +
                theme(axis.text.x = element_text(size=8, angle=45),axis.title.x =element_blank(),
                      axis.title.y=element_blank())+
                facet_wrap(~chemistry,scales = "free")
}
create_plot_10cell_level <- function(x){
        tmp <- ScaleData(tmp,features = x)
        meta <- tmp@meta.data
        expr <- t(tmp@assays$RNA@scale.data)
        df <- do.call(rbind,lapply(unique(meta$sample),function(x){
                sample <- rownames(meta)[meta$sample==x]
                data.frame(sample=x,expr=rollapply(expr[sample,1],10,by=10,mean))
        }))
        df <- data.frame(df)
        df$clinical_outcome <- meta[match(df$sample,meta$sample),'clinical_outcome']
        df$cohort <- meta[match(df$sample,meta$sample),'chemistry']
        df$cohort <- as.factor(df$cohort)
        levels(df$cohort) <- c("BTKi cohort","CAR-T cohort")
        df <- arrange(df,factor(clinical_outcome,levels=c("BTKi-R","Dual-R")))
        lst <- unique(df$sample)
        ggplot(df, aes(factor(sample,levels=lst), expr, fill = clinical_outcome)) +
                geom_boxplot()+
                scale_fill_manual(values=c("peachpuff4","chocolate2"))+
                ylab("Relative expression") +
                theme_bw()+xlab("Clinical outcome") +
                theme(axis.text.x = element_text(size=6, angle=45,h=1),axis.title.x = element_blank())+
                ggtitle(paste0(x," expression"))+
                facet_wrap(~cohort,scales = "free")
}
create_plot_sample_level <- function(x){
        tmp <- ScaleData(tmp,features = x)
        meta <- tmp@meta.data
        expr <- t(tmp@assays$RNA@scale.data)
        expr_mean <- unlist(lapply(unique(meta$sample),function(x){
                mean(expr[rownames(meta)[meta$sample==x],])
        }))
        df <- data.frame(expr=expr_mean,sample=unique(meta$sample))
        df$clinical_outcome <- meta[match(df$sample,meta$sample),'clinical_outcome']
        df$cohort <- meta[match(df$sample,meta$sample),'chemistry']
        df$cohort <- as.factor(df$cohort)
        levels(df$cohort) <- c("Cohort1","Cohort2")
        df <- arrange(df,factor(clinical_outcome,levels=c("BTKi-R","Dual-R")))
        ggplot(df, aes(clinical_outcome, expr, fill = clinical_outcome)) +
                geom_boxplot()+geom_point()+
                scale_fill_manual(values=c("peachpuff4","chocolate2"))+
                ylab("Relative expression") +
                theme_bw()+xlab("Clinical outcome") +
                theme(axis.text.x = element_text(size=6, angle=45,h=1),axis.title.x = element_blank())+
                ggtitle(paste0(x," expression"))+
                #stat_compare_means(comparisons = list(c("BTKi-R","Dual-R")))+
                facet_wrap(~cohort,scales = "free")
}
create_plot_10cell_level('CDK9') ->p1
create_plot_10cell_level('POLR2C') ->p2
p1/p2

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



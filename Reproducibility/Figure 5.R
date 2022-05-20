setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
library(pheatmap)
load("/Users/yanfangfang/BlueYellowColormaps_V1.RData")
b_seu <- readRDS("data/integrated_b_cells.rds")
b_seu$clinical_outcome <- factor(b_seu$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(b_seu$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
subset <- b_seu[,b_seu$clinical_outcome %in% c("BTKi-R","Dual-R")]
good <- names(table(subset$sample))[table(subset$sample)>50]
subset <- subset[,subset$sample %in% good]
DefaultAssay(subset) <- "RNA"
subset$sample[subset$sample=="A3" & subset$chemistry=="5prime"] <- "A3_cohort2"
subset$sample[subset$sample=="A3" & subset$chemistry=="3prime"] <- "A3_cohort1"

plotGene <- function(x){
        subset <- ScaleData(subset,features = x)
        meta <- subset@meta.data
        expr <- t(subset@assays$RNA@scale.data)
        df <- do.call(rbind,lapply(unique(meta$sample),function(x){
                sample <- rownames(meta)[meta$sample==x]
                n <- ceiling(length(sample)/10)
                asplit <- split(sample,cut(1:length(sample),n,label=F))
                do.call(rbind,lapply(asplit,function(y){c(x,round(mean(expr[y,1]),2))}))
        }))
        df <- data.frame(df)
        colnames(df) <- c("sample","expr")
        df$expr <- as.numeric(df$expr)
        df$clinical_outcome <- meta[match(df$sample,meta$sample),'clinical_outcome']
        df <- arrange(df,factor(clinical_outcome,levels=c("BTKi-R","Dual-R")))
        df$cohort <- meta[match(df$sample,meta$sample),'chemistry']
        df$cohort <- as.factor(df$cohort)
        levels(df$cohort) <- c("Cohort1","Cohort2")
        ggplot(df, aes(clinical_outcome, expr, fill = clinical_outcome)) +
                geom_point(size=0.5,aes(color=clinical_outcome),position = position_jitter(seed = 1, width = 0.1))+
                geom_violin(aes(color=clinical_outcome),color="grey",width=0.6, alpha=0.2)+
                scale_color_manual(values=c("peachpuff4","chocolate2"))+
                ylab("Relative expression") +
                theme_bw()+xlab("Clinical outcome") +
                theme(axis.text.x = element_text(size=8, angle=45,h=1),axis.title.x = element_blank())+
                ggtitle(paste0(x," expression"))+
                stat_compare_means(comparisons = list(c("BTKi-R","Dual-R")))+
                facet_wrap(~cohort,scales = "free_y")
}
plotGene("CDK9")

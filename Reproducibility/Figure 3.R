setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
library(pheatmap)
b_seu <- readRDS("data/integrated_b_cells.rds")
good <- names(table(b_seu$sample))[(table(b_seu$sample)>20)]
b_seu <- b_seu[,b_seu$sample %in% good]
b_seu$clinical_outcome <- factor(b_seu$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(b_seu$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")

# Figure 3A ####
DimPlot(b_seu,group.by = 'Phase',cols = c("grey", "orange", "red"))+
        labs(title="Cell cycle phase")
b_seu$tmp <- "Sensitive (BTKi-Fast & BTKi-Slow)"
b_seu$tmp[b_seu$ibrutinib_sensitivity %in% c("R","Dual")] <- "Resistant (BTKi-R & Dual-R)"
DimPlot(b_seu,group.by = 'Phase',split.by="tmp",cols = c("grey", "orange", "red"),ncol=1)+
        NoLegend()+theme(plot.title = element_blank())

# Figure 3B ####
library(dplyr)
df <- b_seu@meta.data %>% group_by(clinical_outcome,patient,sample,Phase) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n))
df <- df[df$Phase!="G1",]
df2 <- df %>% group_by(clinical_outcome,patient,sample) %>% 
        mutate(proliferation=sum(freq))
df2 <- as.data.frame(df2)
df2 <- subset(df2,select=-c(Phase,n,freq))
df2 <- df2[!(duplicated(df2)),]
library(viridis)
library(ggpubr)
my_comparisons <- list(c("BTKi-Fast","BTKi-Slow"),c("BTKi-Fast","BTKi-R"),c("BTKi-Fast","Dual-R"))
ggplot(df2,aes(factor(clinical_outcome,
                      levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")),
               proliferation,label=sample))+
        geom_boxplot()+ylab("Proliferation rate")+
        geom_text(check_overlap = TRUE,
                  position=position_jitter(width=0.15))+
        geom_jitter(aes(color=clinical_outcome), size=2, alpha=0.9)+
        scale_color_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        xlab("Clinical outcome")+
        ggtitle("Proliferation rate, (S+G2M)/Total")+theme_bw()+
        #stat_compare_means(comparisons = my_comparisons)+
        stat_compare_means(label.y=0.85)
# kruskal.test(proliferation ~ clinical_outcomes, data = df2)

# Figure 3C, validate ki67 ####
library(ggplot2)
library(dplyr)
library("readxl")
ki67 <- read_excel("meta/Patient info_scRNA_29 Pt_summary_12302021_vivian.xlsx",sheet="Ki67")
ki67 <- data.frame(t(ki67)[-1,])
colnames(ki67) <- "ki67"
ki67 <- na.omit(ki67)
ki67$ki67 <- as.integer(ki67$ki67)
meta <- read.delim('meta/MCL_meta_0118_2022.txt')
ki67$clinical_outcome <- meta[match(rownames(ki67),meta$patient),'clinical_outcome']
ki67$clinical_outcome <- factor(ki67$clinical_outcome,
                                levels=c("IBN-S","IBN-Slow","IBN-R","Dual"))
levels(ki67$clinical_outcome) <- c("BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
ki67$patient <- rownames(ki67)
library(viridis)
library(ggpubr)
my_comparisons <- list(c("BTKi-Slow","BTKi-R"),c("BTKi-Fast","BTKi-R"))
ggplot(ki67,aes(clinical_outcome,
                ki67,label=patient,fill=clinical_outcome))+
        geom_boxplot()+ylab("Proliferation rate")+
        scale_fill_manual(values=c("steelblue","lightblue2","peachpuff4","chocolate2"))+
        geom_text(check_overlap = TRUE,
                  position=position_jitter(width=0.15))+
        geom_jitter(size=1, alpha=0.6)+
        xlab("Clinical outcome")+
        ggtitle("Proliferation rate (%, Ki67)")+theme_bw()+
        stat_compare_means(comparisons = my_comparisons)+
        stat_compare_means(label.y=130)

setwd("/Users/yanfangfang/Downloads/MW/")
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(dplyr)
# Load data ####
b_seu <- readRDS("data/integrated_b_cells.rds")
b_seu$clinical_outcome <- factor(b_seu$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(b_seu$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
good <- names(table(b_seu$sample))[table(b_seu$sample)>20]
tmp <- b_seu[,b_seu$sample %in% good]
DefaultAssay(tmp) <- "RNA"
tmp$sample[tmp$sample=="A3" & tmp$chemistry=="5prime"] <- "A3_cohort2"
tmp$sample[tmp$sample=="A3" & tmp$chemistry=="3prime"] <- "A3_cohort1"

# HSP90 expression in IBN-R and Dual ####
hsp90 <- c("HSP90AB1","HSP90AA1")
hallmark <- read.delim("/Users/yanfangfang/h.all.v7.4.symbols.gmt",header = F,fill = NA)
rownames(hallmark) <- hallmark$V1
hallmark <- hallmark[,-c(1,2)]
x <- lapply(rownames(hallmark),function(x){
        t <- unname(unlist(hallmark[x,]))
        t[t!=""]})
names(x) <- rownames(hallmark)
myc_v1 <- x[['HALLMARK_MYC_TARGETS_V1']]
myc_v2 <- x[['HALLMARK_MYC_TARGETS_V2']]
myc_v1 <- myc_v1[myc_v1 %in% rownames(tmp)]
myc_v1 <- myc_v1[!(myc_v1 %in% c("HSP90AA1","HSP90AB1"))]
myc_v2 <- myc_v2[myc_v2 %in% rownames(tmp)]
tmp <- ScaleData(tmp,vars.to.regress = 'chemistry',features=c(hsp90,"CDK9",myc_v1,myc_v2))
meta <- tmp@meta.data
plotGene <- function(x){
        df <- data.frame(meta,"expr"=colMeans(tmp@assays$RNA@scale.data[x,,drop=F]))
        df$clinical_outcome <- as.character(df$clinical_outcome)
        df <- df[df$clinical_outcome %in% c("BTKi-R","Dual-R"),]
        my_comparisons <- list(c("BTKi-R","Dual-R"))
        ggplot(df,aes(factor(clinical_outcome,levels=c("BTKi-R","Dual-R")),
                      expr,color=clinical_outcome))+
                geom_point(position = position_jitter(seed = 1, width = 0.2))+
                geom_violin(alpha=0.5)+
                theme_bw()+xlab("Clinical outcome")+ylab("Normalized gene expression")+
                stat_compare_means(comparisons = my_comparisons)+
                scale_color_manual(values=c("peachpuff4","chocolate2"))+
                ggtitle(x)
}
p1 <- plotGene("HSP90AB1")+theme(legend.position = "none")
p2 <- plotGene("HSP90AA1")+theme(legend.position = "none")
# p3 <- plotGene("CDK9")+theme(legend.position = "none")
p4 <- plotGene(myc_v1)+theme(legend.position = "none")+ggtitle("MYC_v1")
p5 <- plotGene(myc_v2)+ggtitle("MYC_v2")
p1|p2|p4|p5


# Figure 4D ####
expr <- t(tmp@assays$RNA@scale.data)
df <- data.frame(tmp@meta.data,"HSP90AB1"=expr[,'HSP90AB1'],"HSP90AA1"=expr[,'HSP90AA1'],
                 "CDK9"=expr[,'CDK9'],"myc_v1"=rowMeans(expr[,myc_v1]),"myc_v2"=rowMeans(expr[,myc_v2]))
asplit <- split(rownames(df), df$sample)
corr_sample <- do.call(rbind, lapply(names(asplit), function(x) {
        cor <- cor(df[asplit[[x]],"myc_v1"], df[asplit[[x]],"HSP90AB1",drop=F])
        c(x,as.numeric(cor),as.character(unique(df[df$sample==x,'clinical_outcome'])))
}))
corr_sample <- data.frame(corr_sample)
colnames(corr_sample) <- c('sample',"corr","clinical_outcome")
corr_sample$corr <- as.numeric(corr_sample$corr)
corr_sample <- corr_sample[corr_sample$clinical_outcome %in% c("BTKi-R","Dual-R"),]
corr_sample$clinical_outcome <- factor(corr_sample$clinical_outcome,levels=c("BTKi-R","Dual-R"))
p1 <- ggplot(corr_sample,aes(clinical_outcome,corr,label=sample))+
        geom_boxplot()+geom_smooth(aes(group=1),color="lightblue")+
        geom_text(check_overlap = TRUE,position=position_jitter(width=0.15))+
        geom_jitter(aes(color=clinical_outcome), size=2, alpha=0.9)+
        scale_color_manual(values=c("peachpuff4","chocolate2"))+
        xlab("Clinical outcome")+ylab("Correlation coefficient")+
        ggtitle("Correlation of HSP90AB1 and MYC_v1 at sample level")+theme_bw()
tmp <- df[df$clinical_outcome %in% c("Dual"),]
tmp$HSP90 <- (tmp$HSP90AB1+tmp$HSP90AA1)/2
p2 <- ggplot(tmp,aes(myc_v1,HSP90,color=sample))+geom_point()+geom_smooth(aes(group=sample))+
        theme_bw()+stat_cor(method = "pearson")+
        ggtitle("Correlation of HSP90AB1 and MYC_v1 at sample level")
p1|p2

# load spliced/unspliced data ####
emat <- readRDS("data/emat.rds")
nmat <- readRDS("data/nmat.rds")
s <- colSums(emat)
u <- colSums(nmat)
p <- data.frame(s,u)
p$pcentage <- apply(p,1,function(x){x[1]/sum(x)})

# Figure 4E ####
library(corrplot)
tmp1 <- df[df$sample=="B1",c("HSP90AB1","HSP90AA1","myc_v1","myc_v2")]
tmp2 <- df[df$sample=="B4",c("HSP90AB1","HSP90AA1","myc_v1","myc_v2")]
corrplot(cor(tmp1), tl.col = "black", tl.srt = 45,col.lim = c(0,1),
         type="lower",diag=FALSE,addCoef.col = "black")
corrplot(cor(tmp2), tl.col = "black", tl.srt = 45,col.lim = c(0,1),
         type="upper",diag=FALSE,addCoef.col = "black")

tmp1 <- df[df$sample=="M0",c("HSP90AB1","HSP90AA1","myc_v1","myc_v2")]
tmp2 <- df[df$sample=="M4",c("HSP90AB1","HSP90AA1","myc_v1","myc_v2")]
corrplot(cor(tmp1), tl.col = "black", tl.srt = 45,col.lim = c(0,1),
         type="lower",diag=FALSE,addCoef.col = "black")
corrplot(cor(tmp2), tl.col = "black", tl.srt = 45,col.lim = c(0,1),
         type="upper",diag=FALSE,addCoef.col = "black")


# HSP90 and co-expressed genes ####
# library(xlsx)
# for (x in c("Dual","IBN-R")){
#         print(x)
#         tmp <- seu[,seu$clinical_outcomes==x]
#         #tmp <- FindVariableFeatures(tmp,nfeatures = 5000)
#         tmp <- ScaleData(tmp,vars.to.regress = 'chemistry',features=rownames(tmp))
#         expr <- tmp@assays$RNA@scale.data
#         meta <- tmp@meta.data
#         good <- names(table(meta$sample))[table(meta$sample)>10]
#         meta <- meta[meta$sample %in% good,]
#         expr <- expr[,rownames(meta)]
#         corr <- do.call(rbind, lapply(rownames(expr), function(x) {
#                 res <- cor.test(expr[x,], expr["HSP90AB1",])
#                 c(res$estimate,res$p.value)
#         }))
#         corr <- data.frame(corr)
#         rownames(corr) <- rownames(expr)
#         colnames(corr) <- c("corr","p_value")
#         corr$adj_p <- p.adjust(corr$p_value)
#         res <- corr[corr$adj_p<0.05,]
#         res <- res[order(res$corr,decreasing = T),]
#         res$gene <- rownames(res)
#         res <- res[,c(4,1:3)]
#         library(xlsx)
#         write.xlsx(res, file="HSP90_coexpressed_genes_significant.xlsx", append=TRUE,
#                    sheetName=paste0("within_",x,"_samples"), row.names=FALSE)
#         corr$gene <- rownames(corr)
#         corr <- na.omit(corr)
#         write.xlsx(corr, file="HSP90_coexpressed_genes_allresults.xlsx", append=TRUE,
#                    sheetName=paste0("within_",x,"_samples"), row.names=FALSE)
# }















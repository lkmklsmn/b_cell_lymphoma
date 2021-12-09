setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
library(pheatmap)
load("/Users/yanfangfang/BlueYellowColormaps_V1.RData")
b_seu <- readRDS("data/integrated_b_cells.rds")
DimPlot(b_seu,group.by = 'clinical_outcomes',
        cols=c("chocolate2","darkgrey","peachpuff4","steelblue","lightblue2"))+
        labs(title="Clinical outcomes") ->p1
DimPlot(b_seu,group.by = 'Phase',cols = c("grey", "orange", "red"))+
        labs(title="Cell cycle phase") ->p2
DimPlot(b_seu,group.by = 'patient',label=T,repel=T)+
        labs(title="Patient")+NoLegend()  ->p3
DimPlot(b_seu,group.by = 'sample',label=T,repel=T)+
        labs(title="Sample")+NoLegend() ->p4
(p1|p2)/(p3|p4)

# Figure 2B, Cell cycle ####
meta <- b_seu@meta.data
library(dplyr)
df <- meta %>% group_by(clinical_outcomes,patient,sample,Phase) %>%
        summarise(n = n()) %>%
        mutate(freq = n / sum(n))
df <- as.data.frame(df)
df$Phase <- factor(df$Phase,levels=c("G2M","S","G1"))
lvls <- unique(df$sample)
# ggplot(df, aes(fill=Phase, y=freq, 
#                x=factor(sample,levels=lvls))) + 
#         facet_wrap(~factor(clinical_outcomes,levels=c("Dual","IBN-R","IBN-Slow","IBN-S","Normal")),scales = "free_x")+
#         geom_bar(position="stack", stat="identity")+
#         scale_fill_manual(values=c("red","orange","grey"))+
#         ggtitle("Cell cycle stage composition")+xlab("Sample")+
#         ylab("fraction (%)")+
#         theme_bw()
df <- df[df$Phase!="G1",]
df2 <- df %>% group_by(clinical_outcomes,patient,sample) %>% 
        mutate(proliferation=sum(freq))
df2 <- as.data.frame(df2)
df2$no.cells <- table(b_seu$sample)[df2$sample]
df2 <- subset(df2,select=-c(Phase,n,freq))
df2 <- df2[!(duplicated(df2)),]
library(viridis)
library(ggpubr)
my_comparisons <- list( c("Dual", "IBN-R"), c("IBN-R", "IBN-Slow"), 
                        c("IBN-Slow", "IBN-S"),c("IBN-S", "Normal") )
ggplot(df2,aes(factor(clinical_outcomes,
                      levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual")),
               proliferation,label=sample))+
        geom_boxplot()+ylab("Proliferation rate")+
        scale_fill_viridis(discrete = TRUE, alpha=0.6) +
        geom_text(check_overlap = TRUE,
                  position=position_jitter(width=0.15))+
        geom_jitter(aes(color=clinical_outcomes), size=2, alpha=0.9)+
        xlab("Therapeutic sensitivity")+
        ggtitle("Proliferation rate, (S+G2M)/Total")+theme_bw()+
        #stat_compare_means(comparisons = my_comparisons)+
        stat_compare_means(label.y=1.3) 
anova <- aov(proliferation ~ no.cells+clinical_outcomes, data = df2)
summary(anova)

# Figure 2C, Outcome-specific DEGs ####
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
de_genes <- c(dual,r,slow,s)

# Heatmap ####
tmp <- b_seu[,b_seu$chemistry=="3prime"]
tmp <- ScaleData(b_seu,features = de_genes)
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
anno_col <- data.frame(sample,"chemistry"=tmp@meta.data[match(sample,tmp$sample),'chemistry'],
                       "clinical_outcomes"=tmp@meta.data[match(sample,tmp$sample),'clinical_outcomes'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col <- arrange(anno_col,factor(clinical_outcomes,levels=c("Dual","IBN-R","IBN-Slow","IBN-S","Normal")))
anno_col$chemistry <- "Cohort1"
annotation_colors = list(
        clinical_outcomes = c(Dual="coral4",'IBN-R'="tomato3",'IBN-Slow'="peachpuff4",
                              'IBN-S'="skyblue4",Normal="grey"),
        chemistry=c(Cohort1="#F8766D")
)
annotation_row <- data.frame("DEGs"=c(rep("Dual-specific",length(dual)),rep("IBN-R-specific",length(r)),
                             rep("IBN-Slow-specific",length(slow)),rep("IBN-S-specific",length(s))))
rownames(annotation_row) <- de_genes
pheatmap(means[,rownames(anno_col)],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         annotation_col = anno_col,show_colnames = F,annotation_colors=annotation_colors,
         show_rownames=F,
         annotation_row = annotation_row,
         main="Outcome-specific DE genes") 

# Figure 2D ####
good <- names(table(b_seu$sample))[table(b_seu$sample)>50]
tmp <- b_seu[,b_seu$sample %in% good]
tmp <- ScaleData(tmp,features=rownames(tmp))
expr <- tmp@assays$RNA@scale.data
meta <- tmp@meta.data
create_plot <- function(gene = "CD52"){
        aframe <- data.frame(meta, expr = expr[gene,])
        aframe <- arrange(aframe,factor(clinical_outcomes,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual")))
        lst <- unique(aframe$sample)
        ggplot(aframe, aes(factor(sample,levels=lst), expr, group = sample, color = clinical_outcomes)) +
                geom_boxplot(outlier.size = 0.2) +
                ylab("Expression levels") +
                theme_bw()+xlab("sample") +
                theme(axis.text.x = element_text(size=8, angle=45))
}
create_plot(dual[3])+ggtitle(paste0("Dual-specific, ",dual[3]))+ theme(legend.position = "none") ->p1
create_plot(r[1])+ggtitle(paste0("IBN-R-specific, ",r[1])) ->p2
create_plot(slow[1])+ggtitle(paste0("IBN-Slow-specific, ",slow[1]))+theme(legend.position = "none") ->p3
create_plot(s[2])+ggtitle(paste0("IBN-S-specific, ",s[2])) ->p4
(p1|p2)/(p3|p4)

# Figure 2E ####
res_dual_r <- read.delim("mixed_model_DEG_Dual_vs_R.txt",sep="\t",row.names = 1)
res_dual_r <- na.omit(res_dual_r)
res_dual_r$adj.p <- p.adjust(res_dual_r$pval_Dual,method="BH")
tmp <- res_dual_r[res_dual_r$adj.p<0.1,]
tmp <- tmp[order(tmp$adj.p,decreasing = F),]
dual_r <-  rownames(tmp)
# heatmap
tmp <- b_seu[,b_seu$clinical_outcomes %in% c("IBN-R","Dual")]
DefaultAssay(tmp) <- "RNA"
tmp$sample[tmp$sample=="A3" & tmp$chemistry=="5prime"] <- "A3_C2"
tmp$sample[tmp$sample=="A3" & tmp$chemistry=="3prime"] <- "A3_C1"
tmp <- ScaleData(tmp,features = rownames(tmp))
asplit_cells <- split(rownames(tmp@meta.data), tmp$sample)
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        s2 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        s3 <- Matrix::rowMeans(tmp@assays$RNA@scale.data[dual_r, sample(unlist(x), 20,replace=T)])
        cbind(s1, s2, s3)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 3)))
anno_col <- data.frame(sample,"sensitivity"=tmp@meta.data[match(sample,tmp$sample),'clinical_outcomes'],
                       "chemistry"=tmp@meta.data[match(sample,tmp$sample),'chemistry'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col$chemistry <- as.factor(anno_col$chemistry)
levels(anno_col$chemistry) <- c("Cohort1","Cohort2")
anno_col$sensitivity <- as.character(anno_col$sensitivity)
anno_col <- arrange(anno_col,factor(sensitivity,levels=c("Dual","IBN-R")))
annotation_colors = list(
        sensitivity = c(Dual="tomato3",'IBN-R'="skyblue4"),
        chemistry=c(Cohort1="#F8766D",Cohort2="#00BFC4")
)
anno_col <- arrange(anno_col,factor(anno_col$chemistry,levels=c("Cohort1","Cohort2")))
pheatmap(means[,rownames(anno_col)],cluster_rows = T, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,gaps_col = 15,
         annotation_col = anno_col,show_colnames = F,show_rownames = F,annotation_colors=annotation_colors,
         main="Dual vs IBN-R in combined cohorts") 

# Figure 2F ####
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
tmp2 <- read.delim("/Users/yanfangfang/Downloads/geneset.txt")
bcr <- tmp2[-1,]

seu <- b_seu
DefaultAssay(seu) <- "RNA"
seu$sample[seu$sample=="A3" & seu$chemistry=="3prime"] <- "A3_cohort1"
seu$sample[seu$sample=="A3" & seu$chemistry=="5prime"] <- "A3_cohort2"

myc_v1 <- myc_v1[myc_v1 %in% rownames(seu)]
myc_v2 <- myc_v2[myc_v2 %in% rownames(seu)]
oxphos <- oxphos[oxphos %in% rownames(seu)]
bcr <- bcr[bcr %in% rownames(seu)]

expr <- seu@assays$RNA@data
m <- mean(as.matrix(expr[oxphos,seu$clinical_outcomes=="Normal"]))
OXPHOS <- colMeans(expr[oxphos,])/m
m <- mean(as.matrix(expr[myc_v1,seu$clinical_outcomes=="Normal"]))
MYC_TARGETS_V1 <- colMeans(expr[myc_v1,])/m
m <- mean(as.matrix(expr[myc_v2,seu$clinical_outcomes=="Normal"]))
MYC_TARGETS_V2 <- colMeans(expr[myc_v2,])/m
m <- mean(as.matrix(expr[bcr,seu$clinical_outcomes=="Normal"]))
KEGG_BCR_Pathway <- colMeans(expr[bcr,])/m
df1 <- data.frame(pathway="OXPHOS",expression=OXPHOS,
                  clinical_outcomes=seu$clinical_outcomes)
df2 <- data.frame(pathway="MYC_TARGETS_V1",expression=MYC_TARGETS_V1,
                  clinical_outcomes=seu$clinical_outcomes)
df3 <- data.frame(pathway="MYC_TARGETS_V2",expression=MYC_TARGETS_V2,
                  clinical_outcomes=seu$clinical_outcomes)
df4 <- data.frame(pathway="KEGG_BCR_Pathway",expression=KEGG_BCR_Pathway,
                  clinical_outcomes=seu$clinical_outcomes)
df <- rbind(df1,df2,df3,df4)
df <- df[df$pathway!="OXPHOS",]
p1 <- ggplot(df, aes(x=factor(clinical_outcomes,levels=c("Dual","IBN-R","IBN-Slow","IBN-S","Normal")), 
               y=expression,fill=pathway)) + 
        geom_boxplot(outlier.shape=NA)+ylim(0,6)+
        xlab("outcome")+ylab("Average pathway score")+
        ggtitle("MYC, OXPHOS, and BCR pathway in different clinical outcomes")+
        theme_bw()+geom_hline(yintercept = 1,linetype="dashed")
p2 <- ggplot(df[df$clinical_outcomes!="Dual",], aes(x=factor(clinical_outcomes,levels=c("Dual","IBN-R","IBN-Slow","IBN-S","Normal")), 
                     y=expression,fill=pathway)) + 
        geom_boxplot(outlier.shape=NA)+ylim(0,6)+
        xlab("outcome")+ylab("Average pathway score")+
        ggtitle("MYC, OXPHOS, and BCR pathway in different clinical outcomes")+
        theme_bw()+geom_hline(yintercept = 1,linetype="dashed")
p1|p2


setwd("/Users/yanfangfang/Downloads/MW/")
library(Seurat)
library(ggplot2)

# Figure S6A ####
infercnv_obj <- readRDS('inferCNV/run.final.infercnv_obj_allcells')
seu <- CreateSeuratObject(infercnv_obj@expr.data)
rm(infercnv_obj)
meta <- read.delim("meta/meta_b_cells.txt")
seu@meta.data <- meta[colnames(seu),]
seu <- NormalizeData(seu)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu,features=rownames(seu))
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.4)
seu <- RunUMAP(seu, dims = 1:10, check_duplicates = FALSE)
seu$clinical_outcome <- factor(seu$clinical_outcomes,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
seu$tmp1 <- "non-A3"
seu$tmp1[seu$sample=="A3_3prime"] <- "A3 in cohort1"
seu$tmp2 <- "non-A3"
seu$tmp2[seu$sample=="A3_5prime"] <- "A3 in cohort2"
p1 <- DimPlot(seu,group.by = 'tmp1',cols=c("#F8766D","lightgrey"))+labs(title="A3 in cohort1")+NoLegend()
p2 <- DimPlot(seu,group.by = 'tmp2',cols=c("#F8766D","lightgrey"))+labs(title="A3 in cohort2")+NoLegend()
p3 <- DimPlot(seu,group.by = 'clinical_outcome',cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Clinical outcome")
p4 <- DimPlot(seu,group.by = 'source')+labs(title="Source")
p1|p2|p3|p4

# Figure S6B ####
require(graphics)
dist.mtx <- as.matrix(dist(seu@reductions$umap@cell.embeddings,method = "euclidean",upper=TRUE))
df <- seu@meta.data
lst <- levels(seu$clinical_outcome)
res <- do.call(rbind,lapply(lst,function(x) {
        c1 <- rownames(df)[df$clinical_outcome=="Normal"]
        c2 <- sample(rownames(df)[df$clinical_outcome==x],length(c1))
        df <- dist.mtx[c1,c2]
        dist <- apply(df,2,mean)
        tmp <- data.frame('dist'=dist,x)
        colnames(tmp) <- c("dist",'clinical_outcome')
        tmp
}))
res <- data.frame(res)
res$sample <- df[match(rownames(res),rownames(df)),'sample']
res$source <- df[match(res$sample,df$sample),'source']
library(ggpubr)
res$clinical_outcome <- factor(res$clinical_outcome,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
ggplot(res[res$source=="PB",], aes(x=clinical_outcome,y=dist,fill=clinical_outcome)) + 
        geom_boxplot(outlier.shape=NA)+
        xlab("clincial outcome")+ylab("Euclidean distance")+
        ggtitle("Genome instability score across different clinical outcome (only PB samples)")+
        theme_bw()+scale_fill_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        stat_compare_means(method = 'anova')+theme(legend.position = 'none') 

# Figure S6C ####
lst <- names(table(df$sample))[table(df$sample)>10]
res <- do.call(rbind,lapply(lst,function(x) {
        c1 <- rownames(df)[df$clinical_outcomes=="Normal"]
        c2 <- rownames(df)[df$sample==x]
        df <- dist.mtx[c1,c2]
        dist <- apply(df,2,mean)
        tmp <- data.frame('dist'=dist,x)
        colnames(tmp) <- c("dist",'sample')
        tmp
}))
res$patient <- df[match(res$sample,df$sample),'patient']
res$source <- df[match(res$sample,df$sample),'source']
res$sample[res$sample=="A3_3prime"] <- "A3_cohort1"
res$sample[res$sample=="A3_5prime"] <- "A3_cohort2"
p1 <- ggplot(res, aes(x=sample, y=dist,fill=patient)) + 
        geom_boxplot(outlier.shape=NA)+
        xlab("outcome")+ylab("Euclidean distance")+
        ggtitle("Genome instability score across samples (colored by patient)")+
        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))+
        theme(legend.position = 'none')
p2 <- ggplot(res, aes(x=sample, y=dist,fill=source)) + 
        geom_boxplot(outlier.shape=NA)+
        xlab("outcome")+ylab("Euclidean distance")+
        ggtitle("Genome instability score across samples (colored by source)")+
        theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust=1))
p1/p2

# Figure S6D ####
b_seu <- readRDS('data/b_seu_umap.rds')
umap <- data.frame(b_seu@reductions$umap@cell.embeddings)
rownames(umap) <- paste0("3prime_",rownames(umap))
umap <- umap[rownames(umap) %in% colnames(seu),]
aframe <- lapply(rownames(umap),function(x) {
  c1 <- colnames(seu)[seu$clinical_outcome=="Normal"]
  df <- dist.mtx[c1,x]
  return(mean(df))
})
umap$score <- unlist(aframe)
p1 <- ggplot(umap,aes(UMAP_1,UMAP_2,color=score))+geom_point()+
  ggtitle("Colored by genome instability score")+
  theme_bw()+
  scale_color_gradient2(low="grey",high="darkblue")

# Figure S6E (spliced ratio vs genome instability score) ####
emat <- readRDS("data/emat.rds")
nmat <- readRDS("data/nmat.rds")
s <- colSums(emat)
u <- colSums(nmat)
p <- data.frame(s,u)
p$pcentage <- apply(p,1,function(x){x[1]/sum(x)})
b_seu$percentage <- p[colnames(b_seu),'pcentage']
df <- b_seu@reductions$umap@cell.embeddings
df <- data.frame(df,b_seu@meta.data)
df <- df[!(is.na(df$percentage)),]
p2 <- ggplot(df,aes(UMAP_1,UMAP_2,color=percentage))+
  geom_point()+theme_bw()+
  scale_color_gradient2(low="grey",mid="lightblue",high="darkblue",midpoint = 0.5)+
  ggtitle("Colored by spliced ratio")

rownames(df) <- paste0("3prime_",rownames(df))
df$score <- umap[rownames(df),'score']
df <- df[!(is.na(df$ibrutinib_sensitivity)),]
df$clinical_outcome <- factor(df$ibrutinib_sensitivity,
                                 levels=c("Normal","S","Slow_responder","R","Dual"))
levels(df$clinical_outcome) <- c("Normal","IBN-S","IBN-Slow","IBN-R","Dual")
p1 <- ggplot(df[df$score<10,],aes(percentage,score,color=clinical_outcome))+geom_point()+
  scale_color_manual(values=c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
  xlab("Spliced ratio")+ylab("Genome instability score")+theme_bw()+
  theme(legend.position = "None")

# Figure S6F ####
lst <- names(table(df$sample))[table(df$sample)>10]
res <- do.call(rbind,lapply(lst,function(x){
  tmp <- df[df$sample==x,]
  t <- cor.test(tmp$score,tmp$percentage)
  return(c(t$p.value,t$estimate))
}))
res <- data.frame(res)
rownames(res) <- lst
colnames(res) <- c("p_value","coefficient")
res$adj_p <- p.adjust(res$p_value,method = "BH")
res$sample <- rownames(res)
res$clinical_outcome <- df[match(res$sample,df$sample),'clinical_outcome']
library(ggrepel)
p2 <- ggplot(res,aes(coefficient,-log10(adj_p),label=sample,color=clinical_outcome))+
  geom_point()+geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  theme_bw()+geom_text_repel(aes(label=sample), size=3)+
  xlab("Correlation coefficient")+scale_y_continuous()+
  scale_color_manual(values=c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))


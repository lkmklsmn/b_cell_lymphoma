setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
seu <- readRDS("data/3prime_b_cells.rds")
seu$clinical_outcome <- factor(seu$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(seu$clinical_outcome) <- c("Normal","IBN-S","IBN-Slow","IBN-R","Dual")
runSeurat <- function(tmp,features){
        tmp <- ScaleData(tmp, features=features)
        tmp <- RunPCA(tmp, features=features)
        tmp <- FindNeighbors(tmp, dims = 1:10)
        tmp <- FindClusters(tmp, resolution = 0.1)
        tmp <- RunUMAP(tmp, dims = 1:10, check_duplicates = FALSE)
        tmp
}

# Figure S3A ####
# Unsupervised highly variable genes
seu <- FindVariableFeatures(seu)
unsuper_genes <- VariableFeatures(seu)
unsuper <- runSeurat(seu,features=unsuper_genes)
p1 <- DimPlot(unsuper,group.by = 'clinical_outcome',
        cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Unsupervised highly variable genes")+NoLegend()

# Supervised Wilcox test
Idents(seu) <- seu$clinical_outcome
stan_df <- FindAllMarkers(seu,only.pos = T)
stan_genes <- unique(stan_df$gene)
stan <- runSeurat(seu,features=stan_genes)
p2 <- DimPlot(stan,group.by = 'clinical_outcome',
              cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Supervised Wilcoxon Rank Sum test (default)")+NoLegend()

# Mixed model
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
mixed_genes <- c(s,slow,r,dual)
mixed <- runSeurat(seu,features=mixed_genes)
p3 <- DimPlot(seu,group.by = 'clinical_outcome',
              cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Supervised Mixed model DEGs")
p1|p2|p3

# Figure S3B ####
good <- names(table(seu$sample))[table(seu$sample)>50]
sub <- seu[,seu$sample %in% good]
tmp1 <- RunPCA(sub, features=s)
p1 <- DimPlot(tmp1,reduction = 'pca',group.by = 'clinical_outcome',
              cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="IBN-S-specific DEGs")+NoLegend()
tmp2 <- RunPCA(sub, features=slow)
p2 <- DimPlot(tmp2,reduction = 'pca',group.by = 'clinical_outcome',
              cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="IBN-Slow-specific DEGs")+NoLegend()
tmp3 <- RunPCA(sub, features=r)
p3 <- DimPlot(tmp3,reduction = 'pca',group.by = 'clinical_outcome',
              cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="IBN-R-specific DEGs")+NoLegend()
tmp4 <- RunPCA(sub, features=dual)
p4 <- DimPlot(tmp4,reduction = 'pca',group.by = 'clinical_outcome',
              cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Dual-specific DEGs")+NoLegend()
(p1|p2)/(p3|p4)

# Figure S3C ####
plot_PCA <- function(tmp){
        df <- data.frame(tmp@reductions$pca@cell.embeddings[,c(1,2)])
        df$sample <- tmp$sample
        df$clinical_outcome <- tmp$clinical_outcome
        df <- arrange(df,factor(clinical_outcome,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual")))
        lst <- unique(df$sample)
        ggplot(df, aes(factor(sample,levels=lst), PC_1, group = sample, fill = clinical_outcome)) +
                geom_boxplot(outlier.size = 0.2) +
                scale_fill_manual(values=c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
                ylab("PC1 coordinates") +
                theme_bw()+xlab("sample") +
                theme(axis.text.x = element_text(size=6,angle=45,hjust=1))+
                theme(legend.position = "none")
}
p1 <- plot_PCA(tmp1)+ggtitle("IBN-S-specific DEGs")+
        geom_hline(yintercept=0, linetype="dashed",color = "steelblue", size=1)
p2 <- plot_PCA(tmp2)+ggtitle("IBN-Slow-specific DEGs")+
        geom_hline(yintercept=0.5, linetype="dashed",color = "lightblue2", size=1)
p3 <- plot_PCA(tmp3)+ggtitle("IBN-R-specific DEGs")+
        geom_hline(yintercept=-0.5, linetype="dashed",color = "peachpuff4", size=1)
p4 <- plot_PCA(tmp4)+ggtitle("Dual-specific DEGs")+
        geom_hline(yintercept=-2.5, linetype="dashed",color = "chocolate2", size=1)
(p1|p2)/(p3|p4)





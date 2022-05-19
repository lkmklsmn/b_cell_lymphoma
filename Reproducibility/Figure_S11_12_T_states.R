setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library("WebGestaltR")
library(pheatmap)
load("/Users/yanfangfang/BlueYellowColormaps_V1.RData")
load("microenvironment/T_cell_states.RData")

# Figure 6A ####
query.projected$chemistry <- "3prime"
query.projected2$chemistry <- "5prime"
merged <- merge(query.projected,y=query.projected2)
merged$sample[merged$sample=="A3" & merged$chemistry=="Cohort1"] <- "A3_c1"
merged$sample[merged$sample=="A3" & merged$chemistry=="Cohort2"] <- "A3_c2"
embed <- rbind(query.projected@reductions$umap@cell.embeddings,query.projected2@reductions$umap@cell.embeddings)
merged[["umap"]] <- CreateDimReducObject(embeddings = embed,key = "UMAP_",assay = DefaultAssay(merged))
stateColors_func <- c("#edbe2a", "#A58AFF", "#00B6EB", "#F8766D", "#53B400", 
                      "#d1cfcc","#FF0000", "#87f6a5", "#e812dd")
names(stateColors_func) <- levels(ref$functional.cluster)
merged$functional.cluster <- factor(merged$functional.cluster,levels=levels(ref$functional.cluster))
p1 <- DimPlot(merged,group.by = 'functional.cluster',label=T,repel=T,cols=stateColors_func)+labs(title="T cell states")+NoLegend()

# Figure 6B ####
library(dplyr)
meta <- merged@meta.data
df <- meta %>% group_by(chemistry,patient,sample,functional.cluster) %>%
        summarise(n = n()) %>% mutate(freq = n / sum(n))
df <- as.data.frame(df)
df$functional.cluster <- factor(df$functional.cluster,levels=levels(ref$functional.cluster))
p2 <- ggplot(df, aes(fill=functional.cluster, y=n, x=sample)) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=stateColors_func)+
        ggtitle("T cell states proportations")+
        xlab("Sample")+ylab("Number of T cells")+theme_bw()+ coord_flip()+
        facet_wrap(~chemistry,scales = "free")

# Figure 6C ####
df <- meta
which.types <- table(df$functional.cluster)> 20
# Calc fold change
df$functional.cluster <- factor(df$functional.cluster)
good <- names(table(df$sample))[table(df$sample)>20]
df <- df[df$sample %in% good,]
lst <- split(df$functional.cluster, df$sample)
freq <- lapply(lst,function(x){
        table(x)/sum(table(x))
})
library(reshape2)
foldchange <- freq$M4/freq$M0
foldchange <- sort(foldchange, decreasing = T)
tb.m <- melt(foldchange)
colnames(tb.m) <- c("Cell_state", "Fold_change")
tb.m <- tb.m[-8,]
ggplot(tb.m, aes(x = Cell_state, y = Fold_change, fill = Cell_state))+
        geom_bar(stat = "identity")+
        scale_fill_manual(values = stateColors_func[tb.m$Cell_state])+ 
        geom_hline(yintercept = 1)+
        scale_y_continuous(trans = "log2") +theme_bw()+ylab("Fold change")+
        theme(axis.text.x = element_text(angle = 45,size=8,hjust=1),axis.title.x=element_blank(),
              legend.position = 'None')+
        ggtitle("T cell states enrichment (M4 vs M0)")

# load data ####
require(dplyr)
require(nlme)
combine <- readRDS("data/combine_w_celltype.rds")
combine <- CellCycleScoring(combine,s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes)

# FMC63-28Z ####
x1 <- paste0("3prime_",colnames(query.projected))
x2 <- paste0("5prime_",colnames(query.projected2))
tmp <- combine[,c(x1,x2)]
df1 <- query.projected@reductions$umap@cell.embeddings
rownames(df1) <- paste0("3prime_",rownames(df1))
df2 <- query.projected2@reductions$umap@cell.embeddings
rownames(df2) <- paste0("5prime_",rownames(df2))
embed <- rbind(df1,df2)
tmp[["umap"]] <- CreateDimReducObject(embeddings = embed,key = "UMAP_",assay = DefaultAssay(tmp))
FeaturePlot(tmp,features = "FMC63-28Z",reduction = 'umap')
x <- combine@assays$RNA@counts['FMC63-28Z',]
summary(x>0)

# CD4 Th1 cells ####
type <- "Th1"
x1 <- paste0("3prime_",colnames(query.projected)[query.projected$functional.cluster==type])
x2 <- paste0("5prime_",colnames(query.projected2)[query.projected2$functional.cluster==type])
tmp <- combine[,c(x1,x2)]
tmp <- tmp[,tmp$ibrutinib_sensitivity %in% c("R","Dual") & tmp$chemistry=="5prime"]
meta <- tmp@meta.data
good <- names(table(meta$sample))[table(meta$sample)>20]
tmp <- tmp[,tmp$sample %in% good]
meta <- meta[meta$sample %in% good,]
meta$ibrutinib_sensitivity <- factor(meta$ibrutinib_sensitivity, levels = c("R", "Dual"))
tmp <- ScaleData(tmp, features = rownames(tmp))
expr <- tmp@assays$RNA@scale.data
sums <- log(rowSums(tmp@assays$RNA@counts))
genes <- names(which(sums > 3))
res_Dual_R <- do.call(rbind, lapply(genes, function(x){
        aframe <- data.frame(meta, expr = expr[x,])
        m1 <- try(lme(expr ~ S.Score + G2M.Score + ibrutinib_sensitivity, random = ~ 1|patient, data = aframe))
        if(class(m1) == "try-error") return(rep(NA, 3))
        as.numeric(coefficients(summary(m1))[4, c(1,4,5)])
}))
rownames(res_Dual_R) <- genes
colnames(res_Dual_R) <- c("coef_Dual", "t_Dual", "pval_Dual")
res_Dual_R <- na.omit(data.frame(res_Dual_R))
res_Dual_R$adj_p <- p.adjust(res_Dual_R$pval_Dual,method="BH")
saveRDS(res_Dual_R,file="microenvironment/res_Dual_R_Th1.rds")
# GSEA 
rank <- res_Dual_R[,'t_Dual',drop=F]
rank <- rank[order(rank$t_Dual,decreasing = T),,drop=F]
rank$gene <- rownames(rank)
rank <- rank[,c(2,1)]
enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                            enrichDatabaseFile = c("/Users/yanfangfang/h.all.v7.4.symbols.gmt"),
                            enrichDatabaseType="genesymbol",
                            interestGene=rank,
                            minNum=5,fdrThr = 0.1,
                            interestGeneType="genesymbol", 
                            outputDirectory="microenvironment/",
                            projectName="Th1_Dual_up_GSEA")
# heatmap
asplit_cells <- split(rownames(tmp@meta.data), tmp$sample)
sig <- rownames(res_Dual_R)[res_Dual_R$coef_Dual>0 & res_Dual_R$adj_p<0.05]
sig <- sig[grep("^RPS|RPL|MT",sig,invert = T)]
gene.use <- sig
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(expr[gene.use, sample(unlist(x), 10)])
        s2 <- Matrix::rowMeans(expr[gene.use, sample(unlist(x), 10)])
        s3 <- Matrix::rowMeans(expr[gene.use, sample(unlist(x), 10)])
        cbind(s1, s2, s3)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 3)))
anno_col <- data.frame(sample,"sensitivity"=tmp@meta.data[match(sample,tmp$sample),'ibrutinib_sensitivity'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col$sensitivity <- as.factor(anno_col$sensitivity)
levels(anno_col$sensitivity) <- c("Dual","IBN_R")
anno_col <- arrange(anno_col,factor(sensitivity,levels=c("Dual","IBN_R")))
annotation_colors = list(sensitivity = c(Dual="chocolate2",IBN_R="peachpuff4"))
pheatmap(means[,rownames(anno_col)],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         annotation_col = anno_col,show_colnames = F,annotation_colors=annotation_colors,
         fontsize_row = 7,
         main="Dual vs IBN-R among Th1 cells")


# CD8 EffectorMemory cells ####
type <- "CD8_EffectorMemory"
x1 <- paste0("3prime_",colnames(query.projected)[query.projected$functional.cluster==type])
x2 <- paste0("5prime_",colnames(query.projected2)[query.projected2$functional.cluster==type])
tmp <- combine[,c(x1,x2)]
tmp <- tmp[,tmp$ibrutinib_sensitivity %in% c("R","Dual") & sub$chemistry=="5prime"]
meta <- tmp@meta.data
good <- names(table(meta$sample))[table(meta$sample)>20]
tmp <- tmp[,tmp$sample %in% good]
meta <- meta[meta$sample %in% good,]
meta$ibrutinib_sensitivity <- factor(meta$ibrutinib_sensitivity, levels = c("R", "Dual"))
tmp <- ScaleData(tmp, features = rownames(tmp))
expr <- tmp@assays$RNA@scale.data
sums <- log(rowSums(tmp@assays$RNA@counts))
genes <- names(which(sums > 3))
res_Dual_R <- do.call(rbind, lapply(genes, function(x){
        aframe <- data.frame(meta, expr = expr[x,])
        m1 <- try(lme(expr ~ S.Score + G2M.Score + ibrutinib_sensitivity, random = ~ 1|patient, data = aframe))
        if(class(m1) == "try-error") return(rep(NA, 3))
        as.numeric(coefficients(summary(m1))[4, c(1,4,5)])
}))
rownames(res_Dual_R) <- genes
colnames(res_Dual_R) <- c("coef_Dual", "t_Dual", "pval_Dual")
res_Dual_R <- na.omit(data.frame(res_Dual_R))
res_Dual_R$adj_p <- p.adjust(res_Dual_R$pval_Dual,method="BH")
saveRDS(res_Dual_R,file="microenvironment/res_Dual_R_CD8_EffectorMemory.rds")
asplit_cells <- split(rownames(tmp@meta.data), tmp$sample)
res_Dual_R <- res_Dual_R[order(res_Dual_R$t_Dual,decreasing = T),]
sig <- rownames(res_Dual_R)[res_Dual_R$coef_Dual>0 & res_Dual_R$adj_p<0.05]
sig <- sig[grep("^RPS|RPL|MT",sig,invert = T)]
gene.use <- sample(sig,80)
means <- do.call(cbind, lapply(asplit_cells, function(x){
        s1 <- Matrix::rowMeans(expr[gene.use, sample(unlist(x), 20)])
        s2 <- Matrix::rowMeans(expr[gene.use, sample(unlist(x), 20)])
        s3 <- Matrix::rowMeans(expr[gene.use, sample(unlist(x), 20)])
        cbind(s1, s2, s3)
}))
sample <- unlist(lapply(names(asplit_cells), function(x) rep(x, 3)))
anno_col <- data.frame(sample,"sensitivity"=tmp@meta.data[match(sample,tmp$sample),'ibrutinib_sensitivity'])
rownames(anno_col) <- colnames(means) <- paste(colnames(means), sample)
anno_col$sensitivity <- as.factor(anno_col$sensitivity)
levels(anno_col$sensitivity) <- c("Dual","IBN_R")
anno_col <- arrange(anno_col,factor(sensitivity,levels=c("Dual","IBN_R")))
annotation_colors = list(sensitivity = c(Dual="chocolate2",IBN_R="peachpuff4"))
pheatmap(means[,rownames(anno_col)],cluster_rows = F, cluster_cols = F, scale = "row",
         breaks = seq(-2, 2, length = length(yellow2blue) + 1), col = yellow2blue,
         annotation_col = anno_col,show_colnames = F,annotation_colors=annotation_colors,
         fontsize_row = 7,
         main="Dual vs IBN-R among CD8_EffectorMemory cells")
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                            enrichDatabase=c("geneontology_Biological_Process"),
                            enrichDatabaseFile = c("/Users/yanfangfang/h.all.v7.4.symbols.gmt"),
                            enrichDatabaseType="genesymbol",
                            interestGene=sig,
                            minNum=5,maxNum = 500,fdrThr = 0.1,
                            interestGeneType="genesymbol",
                            referenceSet ="genome_protein-coding",
                            referenceGeneType="genesymbol", isOutput=TRUE,
                            outputDirectory="microenvironment/",
                            projectName="CD8_EffectorMemory_Dual_up_ORA")
rank <- res_Dual_R[,'t_Dual',drop=F]
rank$gene <- rownames(rank)
rank <- rank[,c(2,1)]
enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                            enrichDatabaseFile = c("/Users/yanfangfang/h.all.v7.4.symbols.gmt"),
                            enrichDatabaseType="genesymbol",
                            interestGene=rank,
                            minNum=5,maxNum = 500,fdrThr = 0.1,
                            interestGeneType="genesymbol", 
                            outputDirectory="microenvironment/",
                            projectName="CD8_EffectorMemory_Dual_up_GSEA")

df <- read.delim("microenvironment/Project_CD8_EffectorMemory_Dual_up_ORA/enrichment_results_CD8_EffectorMemory_Dual_up_ORA.txt")
df <- df[df$FDR<0.01 & df$enrichmentRatio>5,c(2,7,11)]
colnames(df) <- c('geneSet','EnrichmentScore','database')
df <- df[order(df$EnrichmentScore,decreasing = T),]
df <- df[c(1:15),]
df2 <- read.delim("microenvironment/Project_CD8_EffectorMemory_Dual_up_GSEA/enrichment_results_CD8_EffectorMemory_Dual_up_GSEA.txt")
df2$database <- "Hallmark"
df2 <- df2[1:2,c(1,4,12)]
colnames(df2) <- c('geneSet','EnrichmentScore','database')
tmp <- data.frame(rbind(df,df2))
ggplot(tmp,aes(reorder(geneSet, EnrichmentScore),EnrichmentScore))+geom_bar(stat="identity")+
        coord_flip()+facet_wrap(~database,scales = "free_x",ncol=1)+
        theme_bw()+xlab("Enrichment score")+ggtitle("Gene Ontology, Biological Process (GO_BP)")
ggplot(df2,aes(reorder(geneSet, EnrichmentScore),EnrichmentScore))+geom_bar(stat="identity")+
        coord_flip()+
        theme_bw()+xlab("Enrichment score")+ggtitle("Hallmark")

# Myeloid cells ####
myeloid <- readRDS("data/five_myeloid.rds")
myeloid$tmp <- as.character(Idents(myeloid))
myeloid <- myeloid[,myeloid$tmp!="Doublets"]
DimPlot(myeloid,group.by = 'tmp',label=T,repel=T)+labs(title="Myeloid cells")+NoLegend()
# M_MDSC: low in HLA_DRA, high in CD84 & CD33
VlnPlot(myeloid,c("ITGAM","CD14","HLA-DRA","HLA-DRB1"),pt.size=0,ncol = 2)

tmp <- myeloid[,myeloid$tmp=="M-MDSC"]
meta <- tmp@meta.data
good <- names(table(meta$sample))[table(meta$sample)>20]
tmp <- tmp[,tmp$sample %in% good]
meta <- meta[meta$sample %in% good,]
meta$ibrutinib_sensitivity <- factor(meta$ibrutinib_sensitivity, levels = c("R", "Dual"))
tmp <- ScaleData(tmp, features = rownames(tmp))
expr <- tmp@assays$RNA@scale.data
sums <- log(rowSums(tmp@assays$RNA@counts))
genes <- names(which(sums > 3))
res_Dual_R <- do.call(rbind, lapply(genes, function(x){
        aframe <- data.frame(meta, expr = expr[x,])
        m1 <- try(lme(expr ~ S.Score + G2M.Score + ibrutinib_sensitivity, random = ~ 1|patient, data = aframe))
        if(class(m1) == "try-error") return(rep(NA, 3))
        as.numeric(coefficients(summary(m1))[4, c(1,4,5)])
}))
rownames(res_Dual_R) <- genes
colnames(res_Dual_R) <- c("coef_Dual", "t_Dual", "pval_Dual")
res_Dual_R <- na.omit(data.frame(res_Dual_R))
res_Dual_R$adj_p <- p.adjust(res_Dual_R$pval_Dual,method="BH")
saveRDS(res_Dual_R,file="microenvironment/res_Dual_R_m-mdsc.rds")

# Figure S12C #### 
# sig <- rownames(res_Dual_R)[res_Dual_R$coef_Dual>0 & res_Dual_R$adj_p<0.05]
# enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#                             enrichDatabase=c("geneontology_Biological_Process_noRedundant"),
#                             enrichDatabaseFile = c("/Users/yanfangfang/h.all.v7.4.symbols.gmt"),
#                             enrichDatabaseType="genesymbol",
#                             interestGene=sig,
#                             minNum=5,maxNum = 500,fdrThr = 0.05,
#                             interestGeneType="genesymbol",
#                             referenceSet ="genome_protein-coding",
#                             referenceGeneType="genesymbol", isOutput=TRUE,
#                             outputDirectory="microenvironment/",
#                             projectName="M-MDSC_Dual_up_ORA")
df <- read.delim("microenvironment/Project_M_MDSC_Dual_up_ORA/enrichment_results_M_MDSC_Dual_up_ORA.txt")
df <- df[df$FDR<0.05,c('geneSet','description','enrichmentRatio','FDR')]
df <- df[order(df$enrichmentRatio,decreasing = T),]
df <- df[c(grep("*IL|interleukin*",df$geneSet),grep("*IL|interleukin*",df$description)),]
colnames(df) <- c('geneSet','description','EnrichmentScore','FDR')
df[c(1,2),'description'] <- df[c(1,2),'geneSet']
ggplot(df,aes(reorder(description,EnrichmentScore),y=log2(EnrichmentScore))) +
        geom_bar(stat="identity")+ylab("EnrichmentScore")+
        coord_flip()+theme_bw()+xlab("pathway")+
        ggtitle("Dual-R vs BTKi-R in M-MDSC cells")+
        theme(axis.title.y =element_blank())



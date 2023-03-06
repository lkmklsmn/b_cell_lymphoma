# Load R libs ####
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(WebGestaltR)

# Load data ####
data <- read.csv("../data/bulk_counts.csv",row.names=1)
data <- data[rowSums(data)>0,]
group <- read.csv("../data/bulk_group.csv",header = F)[,2]
colnames(data) <- group

# Initiate DESeq2 object ####
condition <- unlist(lapply(colnames(data),function(x){
  strsplit(x,split="_",fixed=T)[[1]][1]}))
condition <- factor(condition,levels=unique(condition))
colData <- data.frame(row.names = colnames(data), condition)
countData <- data[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ condition)
dds <- DESeq(dds)

# Run PCA ####
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
#pdf("PCA.pdf", width=8, height=6)
ggplot(pcaData, aes(PC1, PC2, color=condition,label=name)) +
  geom_point(size=2) + 
  geom_label_repel(size=3,max.overlaps = 15)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot")+theme_bw()
#dev.off()

# Run specific contrasts ####
lapply(levels(condition)[2:9],function(x){
  group1 <- "DMSO"
  group2 <- x
  setwd(path)
  dir.create(paste0(group2,"_vs_",group1))
  setwd(paste0(group2,"_vs_",group1,"/"))
  df <- cbind(data[,paste0(group1,"_",seq(3))],
              data[,paste0(group2,"_",seq(3))])
  res <- results(dds,contrast=c("condition",group2,group1))
  res <- cbind(df[rownames(res),],res)
  res <- na.omit(res[order(res$padj),])
  # write.table(res,file="DEG_results.csv",row.names = T,
  #             col.names = T,quote=F,sep=",")
  
  padj.thres <- 0.05
  logFC.thres <- 0.5
  up <- res[res$padj<padj.thres & res$log2FoldChange>logFC.thres,]
  up <- up[order(up$log2FoldChange,decreasing = T),]
  write.table(up,file="up_DEGs.csv",row.names = T,col.names = T,quote=F,sep=",")
  down <- res[res$padj<padj.thres & res$log2FoldChange< -logFC.thres,]
  down <- down[order(down$log2FoldChange,decreasing = F),]
  #write.table(down,file="down_DEGs.csv",row.names = T,col.names = T,quote=F,sep=",")
  
  # Volcano plot
  res$label <- rownames(res)
  # res$label[!(res$label %in% c(rownames(up)[1:5],rownames(down)[1:5]))] <- NA
  res$label[!(rownames(res) %in% c("CDK9","MYC","HSP90AB1",'DNAJB4'))] <- NA
  highlight_df <- res[res$label %in% c("CDK9","MYC","HSP90AB1",'DNAJB4'),]
  res$sigGene <- abs(res$log2FoldChange) > logFC.thres & res$padj<padj.thres
  ggplot2::ggsave(filename = paste0("../volcano/Volcano_",group2, "_vs_",group1,".pdf"),
                  width=6, height=4,
                  plot = ggplot(res, aes(x = log2FoldChange, y = -log10(padj)))+
                    geom_point(size=2,color='grey') + geom_hline(yintercept = -log10(0.05),linetype='dashed')+
                    geom_vline(xintercept = c(-0.5,0.5),linetype='dashed')+
                    geom_point(data=highlight_df,
                               aes(x=log2FoldChange,y=-log10(padj), color=label),size=2)+
                    geom_label_repel(aes(label=label),size=2)+
                    theme_bw()+ ggtitle(paste0(group2, " vs ",group1)), 
                  device = "pdf")
  # Enrichment analysis
  enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase=c("geneontology_Biological_Process_noRedundant",
                                               "geneontology_Cellular_Component_noRedundant",
                                               "geneontology_Molecular_Function_noRedundant"),
                              enrichDatabaseType="genesymbol",
                              interestGene = rownames(up),
                              minNum=5,maxNum = 500,
                              interestGeneType="genesymbol",fdrThr = 0.05,
                              referenceSet ="genome_protein-coding",
                              referenceGeneType="genesymbol", isOutput=TRUE, 
                              projectName="GO_up")
  enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase=c("geneontology_Biological_Process_noRedundant",
                                               "geneontology_Cellular_Component_noRedundant",
                                               "geneontology_Molecular_Function_noRedundant"),
                              enrichDatabaseType="genesymbol",
                              interestGene = rownames(down),
                              minNum=5,maxNum = 500,
                              interestGeneType="genesymbol",fdrThr = 0.05,
                              referenceSet ="genome_protein-coding",
                              referenceGeneType="genesymbol", isOutput=TRUE, 
                              projectName="GO_down")
  enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabaseFile = 'C:/Users/fyan/Desktop/MCL_RNAseq/h.all.v7.5.1.symbols.gmt',
                              enrichDatabaseType="genesymbol",
                              interestGene = rownames(up),
                              minNum=5,maxNum = 500,
                              interestGeneType="genesymbol",fdrThr = 1,
                              referenceSet ="genome_protein-coding",
                              referenceGeneType="genesymbol", isOutput=TRUE, 
                              projectName="HM_up")
  enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabaseFile = 'C:/Users/fyan/Desktop/MCL_RNAseq/h.all.v7.5.1.symbols.gmt',
                              enrichDatabaseType="genesymbol",
                              interestGene = rownames(down),
                              minNum=5,maxNum = 500,
                              interestGeneType="genesymbol",fdrThr = 1,
                              referenceSet ="genome_protein-coding",
                              referenceGeneType="genesymbol", isOutput=TRUE, 
                              projectName="HM_down")
})

# Boxplot of CDK9
tmp <- data.frame(norm['CDK9',])
colnames(tmp) <- 'expr'
tmp$group <- unlist(lapply(rownames(tmp),function(x){
  strsplit(x,split="_",fixed=T)[[1]][1]}))
tmp$group <- factor(tmp$group,levels=levels(condition))
ggplot(tmp, aes(x=factor(group,levels=levels(condition)), y=expr,color=group)) +
  geom_boxplot()+theme_bw()+ggtitle('CDK9 expression')+xlab("group")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Load gene sets ####
tmp <- read.delim("../data/h.all.v7.5.1.symbols.gmt",header = F,fill = NA)
rownames(tmp) <- tmp$V1
tmp <- tmp[,-c(1,2)]
x <- lapply(rownames(tmp),function(x){
  t <- unname(unlist(tmp[x,]))
  t[t!=""]
})
names(x) <- rownames(tmp)
apoptosis <- x[['HALLMARK_APOPTOSIS']]
myc_v1 <- x[['HALLMARK_MYC_TARGETS_V1']]
myc_v2 <- x[['HALLMARK_MYC_TARGETS_V2']]
nfkb <- x[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']]

genesets <- x[c('HALLMARK_APOPTOSIS',
                'HALLMARK_MYC_TARGETS_V1',
                'HALLMARK_MYC_TARGETS_V2',
                'HALLMARK_TNFA_SIGNALING_VIA_NFKB')]

# Normalize data to DMSO ####
norm <- counts(dds, normalized=TRUE)
tmp_df <- data.frame(norm)
tmp_df <- tmp_df[rowSums(tmp_df) > 100, ]
norm2 <- t(apply(tmp_df,1,function(x){
  log2(x/mean(x[1:3]))
}))  

# Define function for generating a heatmap ####
create_heatmaps <- function(geneset){
  farben <- c(DMSO = "grey",
              `AZD 5nM` = "royalblue1",
              `Zela 0.2uM` = "red4",
              `Tane 0.5uM` = "orangered3",
              `AZD(5)+Tane (0.5)` = "purple4",
              `AZD(5)+Zela (0.2)` = "navy")
  
  tmp_df2 <- get_geneset(genesets[[geneset]])
  anno <- data.frame("group" = as.character(condition))
  rownames(anno) <- colnames(tmp_df)
  
  ok <- which(anno$group %in% c("DMSO", "AZD 5nM", "Zela 0.2uM", "AZD(5)+Zela (0.2)"))
  rowOrd <- order(rowMeans(tmp_df2[, anno$group == 'AZD(5)+Zela (0.2)']))
  
  hm_1 <- pheatmap(tmp_df2[rowOrd, ok],
                   main = geneset,
                   cluster_cols = F, cluster_rows = F,
                   breaks = seq(-1,1,length=100),
                   color = colorRampPalette(c("blue", "black", "yellow"))(100),
                   show_rownames = T, show_colnames = F,
                   border_color = NA,
                   fontsize = 6,
                   annotation_colors = list(group = farben),
                   annotation_col = anno)
  
  ok <- which(anno$group %in% c("DMSO", "AZD 5nM", "Tane 0.5uM", "AZD(5)+Tane (0.5)"))
  rowOrd <- order(rowMeans(tmp_df2[, anno$group == 'AZD(5)+Tane (0.5)']))
  
  hm_2 <- pheatmap(tmp_df2[rowOrd, ok],
                   main = geneset,
                   cluster_cols = F, cluster_rows = F,
                   breaks = seq(-1,1,length=100),
                   color = colorRampPalette(c("blue", "black", "yellow"))(100),
                   show_rownames = T, show_colnames = F,
                   border_color = NA,
                   fontsize = 6,
                   annotation_colors = list(group = farben),
                   annotation_col = anno)
  
  gridExtra::grid.arrange(grobs = list(hm_1[[4]],
                                       hm_2[[4]]), ncol = 2)  
}

# Create heatmap PDFs #### 
lapply(names(genesets), function(geneset){
  p <- create_heatmaps(geneset = geneset)
  ggsave(p,
         filename = paste0("/Users/lukas/Downloads/Bulk_RNAseq_heatmap_",
                           geneset, ".pdf"),
         width = 8, height = 8)
})

# Highlight specific genes ####
plot_gene <- function(gene){
  farben <- c(DMSO = "grey",
              `AZD 2.5nM` = "royalblue1",
              `AZD 5nM` = "royalblue1",
              `Zela 0.2uM` = "red4",
              `Zela 0.4uM` = "red4",
              `Tane 0.5uM` = "orangered3",
              `Tane 1.0uM` = "orangered3",
              `AZD(5)+Tane (0.5)` = "purple4",
              `AZD(5)+Zela (0.2)` = "navy")
  
  aframe <- data.frame(gene = norm2[gene,],
                       "group" = factor(condition,
                                        levels = c("DMSO",
                                                   "AZD 2.5nM",
                                                   "AZD 5nM",
                                                   "Tane 0.5uM",
                                                   "Tane 1.0uM",
                                                   "Zela 0.2uM",
                                                   "Zela 0.4uM",
                                                   "AZD(5)+Tane (0.5)",
                                                   "AZD(5)+Zela (0.2)"))) 
  azd_tane <- c("DMSO",
                "AZD 2.5nM",
                "AZD 5nM",
                "Tane 0.5uM",
                "Tane 1.0uM",
                "AZD(5)+Tane (0.5)")
  azd_zela <- c("DMSO",
                "AZD 2.5nM",
                "AZD 5nM",
                "Zela 0.2uM",
                "Zela 0.4uM",
                "AZD(5)+Zela (0.2)")
  p1 <- ggplot(aframe[aframe$group %in% azd_tane, ],
               aes(group, gene, color = group)) +
    geom_boxplot() + geom_jitter() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(y = "Normalized counts [log2 relative DMSO]",
         title= gene) +
    scale_color_manual(values = farben) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  p2 <- ggplot(aframe[aframe$group %in% azd_zela, ],
               aes(group, gene, color = group)) +
    geom_boxplot() + geom_jitter() +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(y = "Normalized counts [log2 relative DMSO]",
         title= gene) +
    scale_color_manual(values = farben) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  gridExtra::grid.arrange(p1, p2, ncol = 2)  
}

genes <- c("CDK4", "SLC19A1", "PRMT3",
           "BIRC3", "CASP3", "CASP7",
           "MCM2", "HNRNPD", "CDK4", "CDK2", "HSP90AB1",
           "FOSL2", "RELB", "NFKBIA", "CCL4", "PRF1", "MCL1", "BIRC2")
genes <- intersect(genes, rownames(norm2))
lapply(genes, function(gene){
  print(gene)
  p <- plot_gene(gene)
  ggsave(p,
         filename = paste0("/Users/lukas/Downloads/Bulk_RNAseq_example_", gene, ".pdf"))
})

# Show correlations ####
norm <- counts(dds, normalized=TRUE)
tmp_df <- data.frame(norm) + 1
tmp_df <- tmp_df[rowSums(tmp_df) > 100, ]
norm2 <- t(apply(tmp_df,1,function(x){
  log2(x/mean(x[1:3]))
}))

asplit <- split(colnames(norm2), condition)
means <- do.call(cbind, lapply(asplit, function(x) rowMeans(norm2[, x])))
correl <- cor(means, use = "pairwise.complete")
pheatmap(correl, display_numbers = T,
         main = "Gene-wise correlation",
         color = viridisLite::viridis(100))

aframe <- data.frame(means[, c("Tane 1.0uM", "Zela 0.4uM")])
ggplot(aframe, aes(Tane.1.0uM, Zela.0.4uM)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c() +
  labs(title = "Correlation of foldchanges",
       y = "Zela [0.4uM]",
       x = "Tane [1.0uM]") +
  ggpubr::stat_cor() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_minimal()
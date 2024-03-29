path <- "C:/Users/fyan/Desktop/MCL_RNAseq/hisat2"
setwd(path)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(WebGestaltR)
data <- read.csv("../counts.csv",row.names=1)
data <- data[rowSums(data)>0,]
group <- read.csv("../group.csv",header = F)[,2]
colnames(data) <- group

# Deseq2 ####
condition <- unlist(lapply(colnames(data),function(x){
  strsplit(x,split="_",fixed=T)[[1]][1]}))
condition <- factor(condition,levels=unique(condition))
colData <- data.frame(row.names = colnames(data), condition)
countData <- data[, rownames(colData)]
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = colData, 
                              design = ~ condition)
dds <- DESeq(dds)

# PCA ####
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pdf("PCA.pdf", width=8, height=6)
ggplot(pcaData, aes(PC1, PC2, color=condition,label=name)) +
  geom_point(size=2) + 
  geom_label_repel(size=3,max.overlaps = 15)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  ggtitle("PCA plot")+theme_bw()
dev.off()

# specific contrast ####
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
  write.table(res,file="DEG_results.csv",row.names = T,
              col.names = T,quote=F,sep=",")
  
  padj.thres <- 0.05
  logFC.thres <- 0.5
  up <- res[res$padj<padj.thres & res$log2FoldChange>logFC.thres,]
  up <- up[order(up$log2FoldChange,decreasing = T),]
  write.table(up,file="up_DEGs.csv",row.names = T,col.names = T,quote=F,sep=",")
  down <- res[res$padj<padj.thres & res$log2FoldChange< -logFC.thres,]
  down <- down[order(down$log2FoldChange,decreasing = F),]
  write.table(down,file="down_DEGs.csv",row.names = T,col.names = T,quote=F,sep=",")
  
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

# Heatmap of MYC taregts ####
tmp <- read.delim("../h.all.v7.5.1.symbols.gmt",header = F,fill = NA)
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

norm <- counts(dds,normalized=TRUE)
gene_set <- myc_v1
gene_set <- gene_set[gene_set %in% rownames(norm)]
tmp_df <- data.frame(norm[gene_set,])
tmp_df <- tmp_df[rowSums(tmp_df)>100,]
anno <- data.frame("group"=as.character(condition))
rownames(anno) <- colnames(tmp_df)
library(pheatmap)
pheatmap(tmp_df,main='Apoptosis targets',cluster_cols = F,
         breaks = seq(-2,2,length=100),cluster_rows = T,
         show_rownames = F,annotation_col=anno)
# normalize to DMSO
tmp_df2 <- t(apply(tmp_df,1,function(x){
  log2(x/mean(x[1:3]))
}))
anno <- data.frame("group"=as.character(condition))
rownames(anno) <- colnames(tmp_df)
library(pheatmap)
pheatmap(tmp_df2,main='MYC targets',cluster_cols = F,
         breaks = seq(-2,2,length=100),cluster_rows = F,
         show_rownames = F,annotation_col=anno)






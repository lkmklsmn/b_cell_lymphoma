setwd("/Users/yanfangfang/Downloads/MW/")
library(Seurat)
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(monocle3)
library(SeuratWrappers)
library(pheatmap)
library(ggplot2)
library("WebGestaltR")
load("/Users/yanfangfang/BlueYellowColormaps_V1.RData")
load("b_seu_monocle3.RData")

# Figure 6A,B ####
colData(cds)$clinical_outcome <- factor(colData(cds)$ibrutinib_sensitivity,
                                           levels=c("Normal","S","Slow_responder","R","Dual"))
levels(colData(cds)$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
p1 <- plot_cells(cds, color_cells_by = "clinical_outcome", cell_size = 0.5)+ 
        theme(legend.position = "right")+ggtitle("Clinical outcomes")+
        scale_color_manual(values=c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
time_palette <- RColorBrewer::brewer.pal(9,"BuGn")
p2 <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 0.5)+ 
        theme(legend.position = "right")+ggtitle("Pseudotime")+
        scale_colour_gradient(low = time_palette[2], high = time_palette[8])
p1|p2

# Figure 6C ####
# Trajectory 2/4 vs 6/7/8
# cds_sub2 <- choose_graph_segments(cds)
cells <- colnames(cds_sub2)
p1 <- plot_cells(cds[,cells], color_cells_by = "clinical_outcome")+ 
        theme(legend.position = "right")+ggtitle("Clinical outcomes")+
        scale_color_manual(values=c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
gene.use <- rownames(b_seu)[rowSums(b_seu@assays$RNA@counts[,cells])>0.1*length(cells)]
cellweight_sub <- cellWeights[cells,c(2,6)]
rownames(pseudotime) <- rownames(cellWeights)
pseu_sub <- pseudotime[cells,c(2,6)]
# sce2 <- fitGAM(counts = b_seu@assays$RNA@counts[gene.use,cells],
#                pseudotime = pseu_sub,
#                cellWeights = cellweight_sub)
sce2 <- readRDS("sce2.rds")
earlyDERes <- earlyDETest(sce2, knots = c(1, 2))
earlyDERes$adj_pval <- p.adjust(earlyDERes$pvalue,method="BH")
earlyDERes <- earlyDERes[order(earlyDERes$waldStat, decreasing = TRUE),]
earlyDERes <- na.omit(earlyDERes)
sigGene <- rownames(earlyDERes)[earlyDERes$adj_pval<0.05]
# overall plot
cellweight_tmp <- cellWeights[,c(2,6)]
cellweight_tmp <- cellweight_tmp[rowSums(cellweight_tmp)>0,]
pseu_tmp <- pseudotime[rownames(cellweight_tmp),c(2,6)]
# gene.use <- sigGene[1:10]
# gene.use <- c("HLA-DQA2","HLA-DRB5")
gene.use <- sigGene[grep("^HNRNP",sigGene)]
counts_tmp <- b_seu@assays$RNA@counts[gene.use,rownames(cellweight_tmp)]
sce_overall2 <- fitGAM(counts = counts_tmp,pseudotime = pseu_tmp,cellWeights = cellweight_tmp)
tmp <- pseu_tmp
tmp[rownames(cellweight_tmp)[cellweight_tmp[,1]==0],1]<-0
tmp[rownames(cellweight_tmp)[cellweight_tmp[,2]==0],2]<-0
v <- max(tmp[tmp[,1]==tmp[,2],1])
p1 <- plotSmoothers(sce_overall2, counts_tmp, gene = gene.use[1])+
        ggtitle(gene.use[1])+
        geom_vline(xintercept=v, linetype="dashed", size=1,color="steelblue")
p2 <- plotSmoothers(sce_overall2, counts_tmp, gene = gene.use[11])+
        ggtitle(gene.use[11])+
        geom_vline(xintercept=v, linetype="dashed", size=1,color="steelblue")
p1|p2

# Figure 6D (heatmap) ####
cellweight_sub <- cellWeights[cells,c(2,6)]
o2 <- rownames(cellweight_sub)[cellweight_sub[,1]==1 & cellweight_sub[,2]==0]
o6 <- rownames(cellweight_sub)[cellweight_sub[,1]==0 & cellweight_sub[,2]==1]
tmp <- order(pseu_sub[o2,1])
o2 <- names(pseu_sub[o2,1,drop=F][tmp,])
tmp <- order(pseu_sub[o6,2])
o6 <- names(pseu_sub[o6,2,drop=F][tmp,])
n <- 20
expr <- b_seu@assays$RNA@counts[sigGene,c(o2,o6)]
expr1 <- t(apply(expr[,o2],1,function(x){rollapply(x,n,mean,by=n)}))
expr2 <- t(apply(expr[,o6],1,function(x){rollapply(x,n,mean,by=n)}))
expr <- cbind(expr1,expr2)
anno <- data.frame(Trajectory=rep(c("Branch to BTKi-R/Dual-R (traj 2/4)",
                                    "Branch to BTKi-Slow (traj 6/7/8)"),c(ncol(expr1),ncol(expr2))),
                   pseudotime=c(rollapply(pseu_sub[o2,1],n,mean,by=n),
                                rollapply(pseu_sub[o6,2],n,mean,by=n)))
anno$clinical_outcome <- c(rollapply(as.character(colData(cds)[o2,'clinical_outcome']),n,FUN=function(x) x[n],by=n),
                           rollapply(as.character(colData(cds)[o6,'clinical_outcome']),n,FUN=function(x) x[n],by=n))
rownames(anno) <- colnames(expr) <- seq(1:ncol(expr))
out <- pheatmap(expr, cluster_rows = T, cluster_cols = F,scale='row',
                breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
                col = yellow2blue,
                annotation_col = anno,
                show_colnames = F,show_rownames = F)
annotation_colors = list(clinical_outcome = c(Normal="darkgrey",'BTKi-Fast'="steelblue",
                                              'BTKi-Slow'="lightblue2",'BTKi-R'="peachpuff4",
                                              'Dual-R'="chocolate2"))
group <- sort(cutree(out$tree_row, k=2))
annotation_row <- data.frame("group"=c(rep("Down",table(group)[2]),rep("Up",table(group)[1])))
rownames(annotation_row) <- rownames(expr)[out$tree_row$order]
pheatmap(expr[rownames(expr)[out$tree_row$order],], cluster_rows = F, cluster_cols = F,scale='row',
         breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
         col = yellow2blue,
         gaps_col = cumsum(c(ncol(expr1),ncol(expr2))),
         annotation_col = anno[,c("clinical_outcome","pseudotime","Trajectory")],
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         show_colnames = F,show_rownames = F)

# Figure 6E ####
df <- read.delim("trajectory_early_driver_GSEA/Project_ORA_early_driver_2_6_up/enrichment_results_ORA_early_driver_2_6_up.txt")
df <- df[df$database!="pathway_KEGG",]
df <- df[df$FDR<0.01,c('geneSet','description','enrichmentRatio','FDR')]
df <- df[order(df$enrichmentRatio,decreasing = T),]
colnames(df) <- c('geneSet','description','EnrichmentScore','FDR')
df <- df[c(1:6,16:19),]
df[c(1,6),'description'] <- df[c(1,6),'geneSet']
ggplot(df,aes(reorder(description,EnrichmentScore),y=EnrichmentScore)) +
        geom_bar(stat="identity")+
        coord_flip()+theme_bw()+xlab("pathway")+
        ggtitle("Dual-R/BTKi-R (trajectory 2/4) vs BTKi-Slow (trajectory 6/7/8)")+
        theme(axis.title.y =element_blank())

# Figure 6F ####
# Trajectory 2,4
# cds_sub1 <- choose_graph_segments(cds)
cells <- colnames(cds_sub1)
p1 <- plot_cells(cds[,cells], color_cells_by = "clinical_outcome")+ 
        theme(legend.position = "right")+ggtitle("Clinical outcomes")+
        scale_color_manual(values=c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
gene.use <- rownames(b_seu)[rowSums(b_seu@assays$RNA@counts[,cells])>0.1*length(cells)]
cellweight_sub <- cellWeights[cells,c(2,4)]
rownames(pseudotime) <- rownames(cellWeights)
pseu_sub <- pseudotime[cells,c(2,4)]
# sce <- fitGAM(counts = b_seu@assays$RNA@counts[gene.use,cells],
#               pseudotime = pseu_sub,
#               cellWeights = cellweight_sub)
sce <- readRDS("sce.rds")
earlyDERes <- earlyDETest(sce, knots = c(1, 2))
earlyDERes$adj_pval <- p.adjust(earlyDERes$pvalue,method="BH")
earlyDERes <- na.omit(earlyDERes)

tmp <- earlyDERes[earlyDERes$adj_pval<0.05,]
tmp <- tmp[order(tmp$waldStat, decreasing = TRUE),]
sigGene <- rownames(tmp)
plotSmoothers(sce, b_seu@assays$RNA@counts[,cells], gene = sigGene[4])+ggtitle(sigGene[4])
# overall plot
tmp <- cellWeights[,c(2,4)]
tmp <- tmp[rowSums(tmp)!=0,]
gene.use <- sigGene[grep("^HNRNP",sigGene)]
counts_sub <- b_seu@assays$RNA@counts[gene.use,rownames(tmp)]
sce_overall <- fitGAM(counts = counts_sub,
                      pseudotime = pseudotime[rownames(tmp),c(2,4)],
                      cellWeights = cellWeights[rownames(tmp),c(2,4)])
tmp <- pseudotime[rownames(tmp),c(2,4)]
tmp[rownames(cellweight_sub)[cellweight_sub[,1]==0],1]<-0
tmp[rownames(cellweight_sub)[cellweight_sub[,2]==0],2]<-0
v <- max(tmp[tmp[,1]==tmp[,2],1])
p1 <- plotSmoothers(sce_overall2, counts_tmp, gene = gene.use[1])+
        ggtitle(gene.use[1])+
        geom_vline(xintercept=v, linetype="dashed", size=1,color="steelblue")


# Figure 6G ####
sigGene <- rownames(earlyDERes)[earlyDERes$adj_pval<0.000001]
# Pseudotime
common <- rownames(cellweight_sub)[rowSums(cellweight_sub)==2]
traj_2 <- rownames(cellweight_sub) [cellweight_sub[,1]==1]
traj_4 <- rownames(cellweight_sub) [cellweight_sub[,2]==1]
branch_2 <- setdiff(traj_2,common)
branch_4 <- setdiff(traj_4,common)
pseu_sub <- pseudotime[rownames(cellweight_sub),c(2,4)]
o <- common[order(pseu_sub[common,1])]
# only branches
o2 <- branch_2[order(pseu_sub[branch_2,1])]
o4 <- branch_4[order(pseu_sub[branch_4,2])]
o2 <- o2[seq(10,length(o2),10)]
o4 <- o4[seq(10,length(o4),10)]

n <- 10
expr <- b_seu@assays$RNA@counts[sigGene,c(o2,o4)]
expr1 <- t(apply(expr[,o2],1,function(x){rollapply(x,n,mean,by=n)}))
expr2 <- t(apply(expr[,o4],1,function(x){rollapply(x,n,mean,by=n)}))
expr <- cbind(expr1,expr2)
anno <- data.frame(Trajectory=rep(c("Branch to Dual-R (traj2)","Branch to BTKi-R (traj4)"),c(ncol(expr1),ncol(expr2))),
                   pseudotime=c(rollapply(pseu_sub[o2,1],n,mean,by=n),
                                rollapply(pseu_sub[o4,2],n,mean,by=n)))
anno$clinical_outcome <- c(rollapply(as.character(colData(cds)[o2,'clinical_outcome']),n,FUN=function(x) x[n],by=n),
                           rollapply(as.character(colData(cds)[o4,'clinical_outcome']),n,FUN=function(x) x[n],by=n))
rownames(anno) <- colnames(expr) <- seq(1:ncol(expr))
out <- pheatmap(expr, cluster_rows = T, cluster_cols = F,scale='row',
                breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
                col = yellow2blue,
                annotation_col = anno,
                show_colnames = F,show_rownames = F)
annotation_colors = list(clinical_outcome = c(Normal="darkgrey",'BTKi-Fast'="steelblue",
                                              'BTKi-Slow'="lightblue2",'BTKi-R'="peachpuff4",
                                              'Dual-R'="chocolate2"))
group <- sort(cutree(out$tree_row, k=2))
annotation_row <- data.frame("group"=c(rep("Down",table(group)[2]),
                                       rep("Up",table(group)[1])))
rownames(annotation_row) <- rownames(expr)[out$tree_row$order]
pheatmap(expr[rownames(expr)[out$tree_row$order],], cluster_rows = F, cluster_cols = F,scale='row',
         breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
         col = yellow2blue,
         gaps_col = cumsum(c(ncol(expr1),ncol(expr2))),
         annotation_col = anno[,c("clinical_outcome","pseudotime","Trajectory")],
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         show_colnames = F,show_rownames = F)

# Figure 6H ####
# sig_24 <- read.delim("sig_24.txt")
# test <- rownames(sig_24)[sig_24$group=="group2"]
# test <- test[-(grep("^RPL|MT|RPS",test))]
# enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
#                             enrichDatabase=c("geneontology_Biological_Process","pathway_KEGG"),
#                             enrichDatabaseFile = c("/Users/yanfangfang/h.all.v7.4.symbols.gmt"),
#                             enrichDatabaseType="genesymbol",
#                             interestGene=test,
#                             minNum=5,maxNum = 500,
#                             interestGeneType="genesymbol", 
#                             referenceSet ="genome_protein-coding",
#                             referenceGeneType="genesymbol", isOutput=TRUE, 
#                             projectName="ORA_early_driver_2_4_up")
df <- read.delim("trajectory_early_driver_GSEA/Project_ORA_early_driver_2_4_up/enrichment_results_ORA_early_driver_2_4_up.txt")
df$set <- paste0(df$geneSet," ",df$description)
df <- df[df$FDR<0.001,c('geneSet','description','enrichmentRatio','FDR')]
colnames(df) <- c('geneSet','description','EnrichmentScore','FDR')
df <- df[order(df$EnrichmentScore,decreasing = T),]
df <- df[1:10,]
df[4,'description'] <- df[4,'geneSet']
df[8,'description'] <- df[8,'geneSet']
df[9,'description'] <- "RNA splicing via transesterification"
ggplot(df,aes(reorder(description,EnrichmentScore),y=EnrichmentScore)) +
        geom_bar(stat="identity")+
        coord_flip()+theme_bw()+xlab("pathway")+
        ggtitle("Dual-R (trajectory 2) vs BTKi-R (trajectory 4)")+
        theme(axis.title.y =element_blank())




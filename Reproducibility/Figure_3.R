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

# Figure 3A,B ####
colData(cds)$clinical_outcome <- as.factor(colData(cds)$ibrutinib_sensitivity)
levels(colData(cds)$clinical_outcome) <- c("Dual","Normal","IBN-R","IBN-S","IBN-Slow")
p1 <- plot_cells(cds, color_cells_by = "clinical_outcome", cell_size = 1)+ 
        theme(legend.position = "right")+ggtitle("Clinical outcomes")
p2 <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1)+ 
        theme(legend.position = "right")+ggtitle("Pseudotime")
p1|p2

# Figure 3D ####
# Trajectory 2,4 ####
# cds_sub1 <- choose_graph_segments(cds)
cells <- colnames(cds_sub1)
p1 <- plot_cells(cds[,cells], color_cells_by = "clinical_outcome")
gene.use <- rownames(b_seu)[rowSums(b_seu@assays$RNA@counts[,cells])>0.1*length(cells)]
cellweight_sub <- cellWeights[cells,c(2,4)]
rownames(pseudotime) <- rownames(cellWeights)
pseu_sub <- pseudotime[cells,c(2,4)]
# sce <- readRDS("sce.rds")
sce <- fitGAM(counts = b_seu@assays$RNA@counts[gene.use,cells],
              pseudotime = pseu_sub,
              cellWeights = cellweight_sub)
earlyDERes <- earlyDETest(sce, knots = c(1, 2))
earlyDERes$adj_pval <- p.adjust(earlyDERes$pvalue,method="BH")
earlyDERes <- na.omit(earlyDERes)

tmp <- earlyDERes[earlyDERes$adj_pval<=0.05,]
tmp <- tmp[order(tmp$waldStat, decreasing = TRUE),]
sigGene <- rownames(tmp)
plotSmoothers(sce, b_seu@assays$RNA@counts[,cells], gene = sigGene[3])
# overall plot
counts_sub <- b_seu@assays$RNA@counts[sigGene[1:3],rownames(cellweight_sub)]
sce_overall <- fitGAM(counts = counts_sub,pseudotime = pseu_sub,cellWeights = cellweight_sub)
tmp <- pseu_sub
tmp[rownames(cellweight_sub)[cellweight_sub[,1]==0],1]<-0
tmp[rownames(cellweight_sub)[cellweight_sub[,2]==0],2]<-0
v <- max(tmp[tmp[,1]==tmp[,2],1])
p2 <- plotSmoothers(sce_overall, counts_sub, gene = sigGene[3])+ggtitle(sigGene[3])+
        geom_vline(xintercept=v, linetype="dashed", size=1,color="steelblue")
p1|p2

# Figure 3C ####
# DE genes
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
expr <- b_seu@assays$RNA@counts[sigGene,c(o2,o4)]
anno <- data.frame(Trajectory=rep(c("Branch1","Branch2"),
                                  c(length(o2),length(o4))),
                   pseudotime=c(pseu_sub[o2,1],pseu_sub[o4,2]),
                   clinical_outcome=colData(cds)[c(o2,o4),'clinical_outcome'])
rownames(anno) <- colnames(expr)
out <- pheatmap(expr, cluster_rows = T, cluster_cols = F,scale='row',
                breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
                col = yellow2blue,
                gaps_col = cumsum(c(length(o2),length(o4))),
                annotation_col = anno,
                show_colnames = F,show_rownames = F)
group <- sort(cutree(out$tree_row, k=2))
annotation_row <- data.frame("group"=c(rep("group1",table(group)[2]),
                                       rep("group2",table(group)[1])))
rownames(annotation_row) <- rownames(expr)[out$tree_row$order]
pheatmap(expr[rownames(expr)[out$tree_row$order],], cluster_rows = F, cluster_cols = F,scale='row',
         breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
         col = yellow2blue,
         gaps_col = cumsum(c(length(o2),length(o4))),
         annotation_col = anno,
         annotation_row = annotation_row,
         show_colnames = F,show_rownames = F)

# Figure 3E ####
sig_24 <- read.delim("sig_24.txt")
test <- rownames(sig_24)[sig_24$group=="group2"]
test <- test[-(grep("^RPL|MT|RPS",test))]
enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                            enrichDatabase=c("geneontology_Biological_Process","pathway_KEGG"),
                            enrichDatabaseFile = c("/Users/yanfangfang/h.all.v7.4.symbols.gmt"),
                            enrichDatabaseType="genesymbol",
                            interestGene=test,
                            minNum=5,maxNum = 500,
                            interestGeneType="genesymbol", 
                            referenceSet ="genome_protein-coding",
                            referenceGeneType="genesymbol", isOutput=TRUE, 
                            projectName="ORA_early_driver_2_4_up")

# Figure 3F ####
# Trajectory 2,6 ####
# cds_sub2 <- choose_graph_segments(cds)
cells <- colnames(cds_sub2)
p1 <- plot_cells(cds[,cells], color_cells_by = "ibrutinib_sensitivity")
gene.use <- rownames(b_seu)[rowSums(b_seu@assays$RNA@counts[,cells])>0.1*length(cells)]
cellweight_sub <- cellWeights[cells,c(2,6)]
rownames(pseudotime) <- rownames(cellWeights)
pseu_sub <- pseudotime[cells,c(2,6)]
# sce2 <- readRDS("sce2.rds")
sce2 <- fitGAM(counts = b_seu@assays$RNA@counts[gene.use,cells],
               pseudotime = pseu_sub,
               cellWeights = cellweight_sub)
earlyDERes <- earlyDETest(sce2, knots = c(1, 2))
earlyDERes$adj_pval <- p.adjust(earlyDERes$pvalue,method="BH")
earlyDERes <- earlyDERes[order(earlyDERes$waldStat, decreasing = TRUE),]
earlyDERes <- na.omit(earlyDERes)
sigGene <- rownames(earlyDERes)[earlyDERes$adj_pval<0.000001]
# overall plot
cellweight_tmp <- cellWeights[,c(2,6)]
cellweight_tmp <- cellweight_tmp[rowSums(cellweight_tmp)>0,]
pseu_tmp <- pseudotime[rownames(cellweight_tmp),c(2,6)]
counts_tmp <- b_seu@assays$RNA@counts[sigGene[1:3],rownames(cellweight_tmp)]
sce_overall2 <- fitGAM(counts = counts_tmp,pseudotime = pseu_tmp,cellWeights = cellweight_tmp)
tmp <- pseu_tmp
tmp[rownames(cellweight_tmp)[cellweight_tmp[,1]==0],1]<-0
tmp[rownames(cellweight_tmp)[cellweight_tmp[,2]==0],2]<-0
v <- max(tmp[tmp[,1]==tmp[,2],1])
p2 <- plotSmoothers(sce_overall2, counts_tmp, gene = sigGene[1])+ggtitle(sigGene[1])+
        geom_vline(xintercept=v, linetype="dashed", size=1,color="steelblue")
p1|p2

# Heatmap (lineage 2 vs 6) ####
cellweight_sub <- cellWeights[cells,c(2,6)]
o2 <- rownames(cellweight_sub)[cellweight_sub[,1]==1 & cellweight_sub[,2]==0]
o6 <- rownames(cellweight_sub)[cellweight_sub[,1]==0 & cellweight_sub[,2]==1]
tmp <- order(pseu_sub[o2,1])
o2 <- names(pseu_sub[o2,1,drop=F][tmp,])
tmp <- order(pseu_sub[o6,2])
o6 <- names(pseu_sub[o6,2,drop=F][tmp,])
expr <- b_seu@assays$RNA@counts[sigGene,c(o2,o6)]
anno <- data.frame(Trajectory=rep(c("Branch1","Branch2"),c(length(o2),length(o6))),
                   pseudotime=c(pseu_sub[o2,1],pseu_sub[o6,2]),
                   clinical_outcome=colData(cds)[c(o2,o6),'clinical_outcome'])
rownames(anno) <- colnames(expr)
out <- pheatmap(expr, cluster_rows = T, cluster_cols = F,scale='row',
                breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
                col = yellow2blue,
                gaps_col = cumsum(c(length(o2),length(o6))),
                annotation_col = anno,
                show_colnames = F,show_rownames = F)
group <- sort(cutree(out$tree_row, k=2))
annotation_row <- data.frame("group"=c(rep("group1",table(group)[1]),
                                       rep("group2",table(group)[2])))
rownames(annotation_row) <- rownames(expr)[out$tree_row$order]
pheatmap(expr[rownames(expr)[out$tree_row$order],], cluster_rows = F, cluster_cols = F,scale='row',
         breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
         col = yellow2blue,
         gaps_col = cumsum(c(length(o2),length(o6))),
         annotation_col = anno,
         annotation_row = annotation_row,
         show_colnames = F,show_rownames = F)


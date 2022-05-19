# Load R libs ####
library(ggplot2)
library(monocle3)
library(pheatmap)
library(RColorBrewer)
library(Seurat)
library(SingleCellExperiment)
library(tradeSeq)

# Load data ####
load("/Users/lukas/OneDrive/Documents/GitHub/b_cell_lymphoma/outputs/b_seu_monocle3.RData")

# Figure 6A,B ####
colData(cds)$clinical_outcome <- factor(colData(cds)$ibrutinib_sensitivity,
                                        levels=c("Normal","S","Slow_responder","R","Dual"))
levels(colData(cds)$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")

p1 <- plot_cells(cds, color_cells_by = "clinical_outcome", cell_size = 1)+ 
        theme(legend.position = "right") +
        ggtitle("Clinical outcomes") +
        scale_color_manual(values = c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))

time_palette <- RColorBrewer::brewer.pal(9,"BuGn")
p2 <- plot_cells(cds, color_cells_by = "pseudotime", cell_size = 1) + 
        theme(legend.position = "right") +
        ggtitle("Pseudotime") +
        scale_colour_gradient(low = time_palette[2], high = time_palette[8])
p1|p2

# Figure 3D ####
# Trajectory 2,4
# cds_sub1 <- choose_graph_segments(cds)
cells <- colnames(cds_sub1)
p1 <- plot_cells(cds[,cells], color_cells_by = "clinical_outcome") + 
        theme(legend.position = "right") +
        ggtitle("Clinical outcomes") +
        scale_color_manual(values = c("lightgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
gene.use <- rownames(b_seu)[rowSums(b_seu@assays$RNA@counts[,cells]) > 0.1*length(cells)]
cellweight_sub <- cellWeights[cells,c(2,4)]
rownames(pseudotime) <- rownames(cellWeights)
pseu_sub <- pseudotime[cells,c(2,4)]
sce <- fitGAM(counts = b_seu@assays$RNA@counts[gene.use, cells],
              pseudotime = pseu_sub,
              cellWeights = cellweight_sub)
#sce <- readRDS("/Users/lukas/OneDrive/Documents/GitHub/b_cell_lymphoma/outputs/sce.rds")
earlyDERes <- earlyDETest(sce, knots = c(1, 2))
#earlyDERes <- readRDS("/Users/lukas/OneDrive/Documents/GitHub/b_cell_lymphoma/outputs/earlyDERes.rds")
earlyDERes$adj_pval <- p.adjust(earlyDERes$pvalue, method = "BH")
earlyDERes <- na.omit(earlyDERes)

tmp <- earlyDERes[earlyDERes$adj_pval<0.05,]
tmp <- tmp[order(tmp$waldStat, decreasing = TRUE),]
sigGene <- rownames(tmp)
plotSmoothers(sce, b_seu@assays$RNA@counts[,cells], gene = "HSP90AB1")

plot_cells(cds[, rownames(sce@colData)], color_cells_by = "pseudotime", reduction_method = "UMAP")

tmp <- monocle3::choose_graph_segments(cds = cds, starting_pr_node = 1, ending_pr_nodes = c(1,2))
DimPlot(b_seu[,colnames(tmp)])

cells <- choose_graph_segments(cds)
aframe <- data.frame(gene = b_seu@assays$RNA@data["HSP90AB1", cells],
                     cds[, cells]@colData,
                     sce$slingshot)
aframe$pt <- aframe$pseudotime.V1
aframe$pt[aframe$cellWeights.Y_4 == 0] <- (-1)*aframe$pt[aframe$cellWeights.Y_4 == 0]

aframe$pt_rank <- NA
aframe$pt_rank[aframe$cellWeights.Y_4 == 0] <- -rank((-1)*aframe$pt[aframe$cellWeights.Y_4 == 0])
aframe$pt_rank[aframe$cellWeights.Y_4 == 1] <- rank(aframe$pt[aframe$cellWeights.Y_4 == 1])

ggplot(aframe, aes(pt_rank, gene, color = pt_rank)) +
        facet_wrap(~ cellWeights.Y_4 == 1, scales = "free_x") +
        geom_point() + geom_smooth() +
        scale_color_gradient2(low = "red", mid = "black", high = "blue", midpoint = 0) +
        theme_bw()

tmp <- aframe[aframe$cellWeights.Y_4 == 1,]
tmp <- tmp[order(tmp$pt),]
tmp$X <- tmp$pt_rank/nrow(tmp)

library(mgcv)

afit <- gam(formula = gene ~ s(X, k = 5), data = tmp)
preds <- predict(afit, data.frame(X = seq(0, 1, length = 20)))
plot(preds)

cells <- rownames(sce$slingshot)[sce$slingshot$cellWeights.Y_4 == 1]
pt <- 

get_pt <- function(cells){
        sce$slingshot[cells, "pseudotime.V1"]
}
get_fitted_vals <- function(gene = "HSP90AB1", cells = cells, pt = pt, n_pts = 20){
        expr <- b_seu@assays$RNA@data[gene, cells]
        aframe <- data.frame(expr, pt)
        aframe$pt_rank <- rank(aframe$pt)
        aframe$pt_rank <- aframe$pt_rank/nrow(aframe)
        afit <- gam(formula = expr ~ s(pt_rank, k = 5), data = aframe)
        preds <- predict(afit, data.frame(pt_rank = seq(0, 1, length = n_pts)))
        return(preds)
}

resis <- get_fitted_vals(gene = "HSP90AB1",
                         cells = rownames(sce$slingshot)[sce$slingshot$cellWeights.Y_4 == 1],
                         pt = get_pt(rownames(sce$slingshot)[sce$slingshot$cellWeights.Y_4 == 1]))

sens <- get_fitted_vals(gene = "HSP90AB1",
                         cells = rownames(sce$slingshot)[sce$slingshot$cellWeights.Y_4 == 0],
                         pt = get_pt(rownames(sce$slingshot)[sce$slingshot$cellWeights.Y_4 == 0]))

aframe <- data.frame(resis, sens, pt = 1:length(resis))
aframe <- reshape2::melt(aframe, id.vars = "pt")
ggplot(aframe, aes(pt, value, group = variable, color = variable)) +
        geom_smooth() +
        theme_classic()

# overall plot
tmp <- cellWeights[,c(2,4)]
tmp <- tmp[rowSums(tmp)!=0,]
counts_sub <- b_seu@assays$RNA@counts[sigGene[1:3],rownames(tmp)]
sce_overall <- fitGAM(counts = counts_sub,
                      pseudotime = pseudotime[rownames(tmp),c(2,4)],
                      cellWeights = cellWeights[rownames(tmp),c(2,4)])
tmp <- pseudotime[rownames(tmp),c(2,4)]
tmp[rownames(cellweight_sub)[cellweight_sub[,1]==0],1]<-0
tmp[rownames(cellweight_sub)[cellweight_sub[,2]==0],2]<-0
v <- max(tmp[tmp[,1]==tmp[,2],1])
p2 <- plotSmoothers(sce_overall, counts_sub, gene = sigGene[3])+
        ggtitle(sigGene[3])+
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
annotation_colors = list(
        clinical_outcome = c(Normal="darkgrey",'IBN-S'="steelblue",
                              'IBN-Slow'="lightblue2",'IBN-R'="peachpuff4",
                              Dual="chocolate2")
)
pheatmap(expr[rownames(expr)[out$tree_row$order],], cluster_rows = F, cluster_cols = F,scale='row',
         breaks = seq(-1, 1, length = length(yellow2blue) + 1), 
         col = yellow2blue,
         gaps_col = cumsum(c(length(o2),length(o4))),
         annotation_col = anno,
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         show_colnames = F,show_rownames = F)

# Figure 3E ####
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


# Figure 3F ####
# Trajectory 2,6 ####
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
counts_tmp <- b_seu@assays$RNA@counts[sigGene[1:10],rownames(cellweight_tmp)]
sce_overall2 <- fitGAM(counts = counts_tmp,pseudotime = pseu_tmp,cellWeights = cellweight_tmp)
tmp <- pseu_tmp
tmp[rownames(cellweight_tmp)[cellweight_tmp[,1]==0],1]<-0
tmp[rownames(cellweight_tmp)[cellweight_tmp[,2]==0],2]<-0
v <- max(tmp[tmp[,1]==tmp[,2],1])
p2 <- plotSmoothers(sce_overall2, counts_tmp, gene = sigGene[8])+ggtitle(sigGene[8])+
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
         annotation_colors = annotation_colors,
         show_colnames = F,show_rownames = F)


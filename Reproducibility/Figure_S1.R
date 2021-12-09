setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
b_seu <- readRDS("data/integrated_b_cells.rds")

# Figure S1A ####
good <- names(table(b_seu$sample))[table(b_seu$sample)>10]
b_seu <- b_seu[,b_seu$sample %in% good]
b_seu$sample[b_seu$sample == "A3" & b_seu$chemistry == "3prime"] <- "A3 in Cohort1"
b_seu$sample[b_seu$sample == "A3" & b_seu$chemistry == "5prime"] <- "A3 in Cohort2"
meta <- b_seu@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1","Normal_2"))] <- "no"
b_seu <- NormalizeData(b_seu)
b_seu <- FindVariableFeatures(b_seu)
b_seu <- ScaleData(b_seu, features = VariableFeatures(b_seu))
expr <- b_seu@assays$RNA@scale.data
asplit <- split(rownames(meta), meta$sample)
pseudobulks <- do.call(cbind, lapply(asplit, function(x) rowMeans(expr[,x])))
num_of_cells <- unlist(lapply(asplit, length))
pca <- prcomp(t(pseudobulks))
udata <- umap::umap(pca$x, n_neighbors = 5)
num_of_cells <- num_of_cells[colnames(pseudobulks)]
meta_pseudo <- meta[match(colnames(pseudobulks), meta$sample),]
aframe <- data.frame(udata$layout,num_of_cells,meta_pseudo, pca$x)
colnames(aframe)[c(1,2)] <- c("UMAP1","UMAP2")
aframe$A3 <- NA
aframe$A3[aframe$sample %in% c("A3 in Cohort1","A3 in Cohort2")] <- aframe$sample
aframe$Cohort <- "Cohort1"
aframe$Cohort[aframe$chemistry=="5prime"] <- "Cohort2"
p1 <- ggplot(aframe, aes(UMAP1, UMAP2, color = Cohort)) +
        geom_point() + theme_bw() +ggtitle("Cohort")+theme(legend.position = c(0.3, 0.5))
p2 <- ggplot(aframe, aes(UMAP1, UMAP2, color = A3)) +
        geom_point() + theme_bw() +ggtitle("A3 sample")+theme(legend.position = c(0.3, 0.5))
p3 <- ggplot(aframe, aes(UMAP1, UMAP2, color = clinical_outcomes)) +
        geom_point() + theme_bw() +ggtitle("clinical outcomes")+theme(legend.position = c(0.3, 0.5))
p4 <- ggplot(aframe, aes(UMAP1, UMAP2, color = tumor)) +
        geom_point() + theme_bw()+ggtitle("tumor")+theme(legend.position = c(0.3, 0.5))
p5 <- ggplot(aframe, aes(UMAP1, UMAP2, color = log(nCount_RNA))) +
        geom_point() + theme_bw()+ggtitle("log(nCount_RNA)")+theme(legend.position = c(0.3, 0.5))
p6 <- ggplot(aframe, aes(UMAP1, UMAP2, color = log(num_of_cells))) +
        geom_point() + theme_bw()+ggtitle("log(number of cells)")+theme(legend.position = c(0.3, 0.5))
p7 <- ggplot(aframe, aes(UMAP1, UMAP2, color = source)) +
        geom_point() + theme_bw()+ggtitle("source")+theme(legend.position = c(0.3, 0.5))
(p1|p2|p7)/(p3|p4|p5)

# Figure S2B ####
combined <- readRDS(file="data/combine_w_celltype.rds")
combined <- combined[,combined$celltype!="Eryth"]
combined$Cohort <- "Cohort1"
combined$Cohort[combined$chemistry=="5prime"] <- "Cohort2"
combined$tmp1 <- NA
combined$tmp1[combined$sample == "A3" & combined$chemistry == "3prime"] <- "A3 in Cohort1"
combined$tmp2 <- NA
combined$tmp2[combined$sample == "A3" & combined$chemistry == "5prime"] <- "A3 in Cohort2"
p1 <- DimPlot(combined,group.by = 'Cohort')+labs(title="Cohort")
p2 <- DimPlot(combined,group.by = 'celltype',cols=c("chocolate2","yellowgreen","wheat3","cyan3",
                                                    "khaki","pink","peachpuff4"))+labs(title="Inferred cell type")
p3 <- DimPlot(combined,group.by = 'tmp1')+labs(title="A3 in Cohort1")
p4 <- DimPlot(combined,group.by = 'tmp2')+labs(title="A3 in Cohort2")
p5 <- FeaturePlot(combined, c("CD8B","IL7R","LTB","GNLY"),ncol=2)
((p1|p2)/(p3|p4))|p5

integrated <- readRDS(file="data/integrated.rds")
DefaultAssay(integrated) <- "RNA"
integrated <- integrated[,integrated$celltype!="Eryth"]
integrated$Cohort <- "Cohort1"
integrated$Cohort[integrated$chemistry=="5prime"] <- "Cohort2"
integrated$tmp1 <- NA
integrated$tmp1[integrated$sample == "A3" & integrated$chemistry == "3prime"] <- "A3 in Cohort1"
integrated@meta.data[sample(colnames(integrated)[integrated$tmp2=="A3 in Cohort2"],300),"tmp1"] <- "A3 in Cohort1"
integrated$tmp2 <- NA
integrated$tmp2[integrated$sample == "A3" & integrated$chemistry == "5prime"] <- "A3 in Cohort2"
p1 <- DimPlot(integrated,group.by = 'Cohort')+labs(title="Cohort")
p2 <- DimPlot(integrated,group.by = 'celltype',cols=c("chocolate2","yellowgreen","wheat3","cyan3",
                                                    "khaki","pink","peachpuff4"))+labs(title="Inferred cell type")
p3 <- DimPlot(integrated,group.by = 'tmp1')+labs(title="A3 in Cohort1")
p4 <- DimPlot(integrated,group.by = 'tmp2')+labs(title="A3 in Cohort2")
p5 <- FeaturePlot(integrated, c("CD8B","IL7R","LTB","GNLY"),ncol=2)
((p1|p2)/(p3|p4))|p5

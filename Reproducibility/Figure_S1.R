setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
combined <- readRDS("data/combine_w_celltype.rds")
combined <- combined[,combined$celltype!="Eryth"]
combined$Cohort <- "Cohort1"
combined$Cohort[combined$chemistry=="5prime"] <- "Cohort2"
combined$tmp1 <- "non-A3"
combined$tmp1[combined$sample == "A3" & combined$chemistry == "Cohort1"] <- "A3 in Cohort1"
combined$tmp2 <- "non-A3"
combined$tmp2[combined$sample == "A3" & combined$chemistry == "Cohort2"] <- "A3 in Cohort2"
combined$clinical_outcome <- factor(combined$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(combined$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")

# combined <- combined[,combined$celltype=="Tumor B"]
# Figure S1A ####
good <- names(table(combined$sample))[table(combined$sample)>10]
combined <- combined[,combined$sample %in% good]
combined$sample[combined$sample == "A3" & combined$chemistry == "3prime"] <- "A3 in Cohort1"
combined$sample[combined$sample == "A3" & combined$chemistry == "5prime"] <- "A3 in Cohort2"
meta <- combined@meta.data
meta$tumor <- "yes"
meta$tumor[which(meta$patient %in% c("Normal_1","Normal_2"))] <- "no"
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined, features = VariableFeatures(combined))
expr <- combined@assays$RNA@scale.data
asplit <- split(rownames(meta), meta$sample)
pseudobulks <- do.call(cbind, lapply(asplit, function(x) rowMeans(expr[,x])))
num_of_cells <- unlist(lapply(asplit, length))
pca <- prcomp(t(pseudobulks))
udata <- umap::umap(pca$x, n_neighbors = 5)
num_of_cells <- num_of_cells[colnames(pseudobulks)]
meta_pseudo <- meta[match(colnames(pseudobulks), meta$sample),]
aframe <- data.frame(udata$layout,num_of_cells,meta_pseudo, pca$x)
colnames(aframe)[c(1,2)] <- c("UMAP1","UMAP2")
aframe$A3 <- "non-A3"
aframe$A3[aframe$sample=="A3 in Cohort1"] <- "A3 in Cohort1"
aframe$A3[aframe$sample=="A3 in Cohort2"] <- "A3 in Cohort2"
aframe$Cohort <- "Cohort1"
aframe$Cohort[aframe$chemistry=="5prime"] <- "Cohort2"
p1 <- ggplot(aframe, aes(UMAP1, UMAP2, color = Cohort)) +
        geom_point() + theme_bw() +ggtitle("Cohort")
p2 <- ggplot(aframe, aes(UMAP1, UMAP2, color = A3)) +
        geom_point() + theme_bw() +ggtitle("A3 sample")+
        scale_color_manual(values=c("#F8766D","#00BFC4","lightgrey"))
p3 <- ggplot(aframe, aes(UMAP1, UMAP2, color = clinical_outcome)) +
        geom_point() + theme_bw() +ggtitle("clinical outcome")+
        scale_color_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
p4 <- ggplot(aframe, aes(UMAP1, UMAP2, color = tumor)) +
        geom_point() + theme_bw()+ggtitle("tumor")
p5 <- ggplot(aframe, aes(UMAP1, UMAP2, color = log(nCount_RNA))) +
        geom_point() + theme_bw()+ggtitle("log(nCount_RNA)")
p6 <- ggplot(aframe, aes(UMAP1, UMAP2, color = log(num_of_cells))) +
        geom_point() + theme_bw()+ggtitle("log(number of cells)")
p7 <- ggplot(aframe, aes(UMAP1, UMAP2, color = source)) +
        geom_point() + theme_bw()+ggtitle("source")
(p1/p3)|(p2/p4)|(p5/p6)

# Figure S1BC ####
p1 <- DimPlot(combined,group.by = 'Cohort')+labs(title="Cohort")
p2 <- DimPlot(combined,group.by = 'celltype',cols=c("chocolate3","yellowgreen","wheat3","cyan3",
                                                    "khaki","pink","peachpuff4"))+labs(title="Inferred cell type")
combined$tmp1 <- "non-A3"
combined$tmp1[combined$sample=="A3 in Cohort1"] <- "A3 in Cohort1"
combined$tmp2 <- "non-A3"
combined$tmp2[combined$sample=="A3 in Cohort2"] <- "A3 in Cohort2"
p3 <- DimPlot(combined,group.by = 'tmp1',cols=c("steelblue",'lightgrey'))+labs(title="A3 in Cohort1")
p4 <- DimPlot(combined,group.by = 'tmp2',cols=c("steelblue",'lightgrey'))+labs(title="A3 in Cohort2")
p5 <- FeaturePlot(combined, c("CD8B","IL7R"),ncol=1)
((p1/p3)|(p2/p4))|p5

# Figure S1DE ####
integrated <- readRDS(file="data/integrated.rds")
DefaultAssay(integrated) <- "RNA"
integrated <- integrated[,integrated$celltype!="Eryth"]
integrated$Cohort <- "Cohort1"
integrated$Cohort[integrated$chemistry=="5prime"] <- "Cohort2"
integrated$tmp1 <- "non-A3"
integrated$tmp1[integrated$sample == "A3" & integrated$Cohort == "Cohort1"] <- "A3 in Cohort1"
integrated$tmp2 <- "non-A3"
integrated$tmp2[integrated$sample == "A3" & integrated$Cohort == "Cohort2"] <- "A3 in Cohort2"
integrated@meta.data[sample(colnames(integrated)[integrated$tmp2=="A3 in Cohort2"],300),"tmp1"] <- "A3 in Cohort1"
p1 <- DimPlot(integrated,group.by = 'Cohort')+labs(title="Cohort")
p2 <- DimPlot(integrated,group.by = 'celltype',cols=c("chocolate2","yellowgreen","wheat3","cyan3",
                                                    "khaki","pink","peachpuff4"))+labs(title="Inferred cell type")
p3 <- DimPlot(integrated,group.by = 'tmp1',cols=c("steelblue",'lightgrey'))+labs(title="A3 in Cohort1")
p4 <- DimPlot(integrated,group.by = 'tmp2',cols=c("steelblue",'lightgrey'))+labs(title="A3 in Cohort2")
p5 <- FeaturePlot(integrated, c("CD8B","IL7R"),ncol=1)
((p1/p3)|(p2/p4))|p5

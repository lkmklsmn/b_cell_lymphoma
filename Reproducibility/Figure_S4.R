setwd("/Users/yanfangfang/Downloads/MW/")
library(ggplot2)

# Figure 2E, chr12,chr17 mutated genes ####
infer <- readRDS("inferCNV/run.final.infercnv_obj_allcells")
expr <- infer@expr.data
# Check chr 12 ####
geneinfo <- (infer@gene_order)
ok <- which(geneinfo$chr == "chr12")
aframe <- data.frame(expr[ok,])
means <- do.call(cbind, lapply(infer@observation_grouped_cell_indices, function(x) rowMeans(expr[ok, x])))
means <- data.frame(means)
means$position <- geneinfo$start[ok]
tmp <- reshape2::melt(means, id.vars = "position")
ggplot(tmp, aes(position, value, group = variable, color = variable)) +
        geom_point() +
        geom_hline(yintercept = 1) +
        geom_vline(xintercept = 7500000, color = "red") +
        theme_classic()
hits <- rownames(geneinfo)[geneinfo$start < 7500000 & geneinfo$chr == "chr12"]
colnames(tmp)[2] <- "sample"
meta <- read.delim("meta/MCL_meta_0118_2022.txt")
meta$sample[meta$sample=="A3_c1"] <- "A3_3prime"
meta$sample[meta$sample=="A3_c2"] <- "A3_5prime"
tmp$clinical_outcome <- meta[match(tmp$sample,meta$sample),'clinical_outcome']
tmp$clinical_outcome <- factor(tmp$clinical_outcome,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
levels(tmp$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
ggplot(tmp[tmp$position < 30000000,], aes(position, value, group = sample, color = clinical_outcome)) +
        geom_point() + ylab("Copy number variation")+
        geom_hline(yintercept = 1) +
        theme_bw()+ggtitle("Chromosome 12")+
        scale_color_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
tmp2 <- geneinfo[rownames(geneinfo) %in% hits,]
tmp2$gene <- rownames(tmp2)
library(ggrepel)
ggplot(tmp2, aes(x = start, xend = stop, y=gene, yend = gene)) +
        geom_segment(linetype = 'solid')+xlab("chrosomose position")+theme_bw()+
        geom_label_repel(aes(label=gene))

# Check chr 17 ####
ok <- which(geneinfo$chr == "chr17")
aframe <- data.frame(expr[ok,])
means <- do.call(cbind, lapply(infer@observation_grouped_cell_indices, function(x) rowMeans(expr[ok, x])))
means <- data.frame(means)
means$position <- geneinfo$start[ok]
tmp <- reshape2::melt(means, id.vars = "position")
ggplot(tmp, aes(position, value, group = variable, color = variable)) +
        geom_point() +
        geom_hline(yintercept = 1) +
        geom_vline(xintercept = 7500000, color = "red") +
        theme_classic()
hits <- rownames(geneinfo)[geneinfo$start < 7500000 & geneinfo$chr == "chr17"]
colnames(tmp)[2] <- "sample"
meta <- read.delim("meta/MCL_meta_0118_2022.txt")
meta$sample[meta$sample=="A3_c1"] <- "A3_3prime"
meta$sample[meta$sample=="A3_c2"] <- "A3_5prime"
tmp$clinical_outcome <- meta[match(tmp$sample,meta$sample),'clinical_outcome']
tmp$clinical_outcome <- factor(tmp$clinical_outcome,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
levels(tmp$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
ggplot(tmp[tmp$position < 30000000,], aes(position, value, group = sample, color = clinical_outcome)) +
        geom_point() + ylab("Copy number variation")+
        geom_hline(yintercept = 1) +
        theme_bw()+ggtitle("Chromosome 17")+
        scale_color_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))
tmp2 <- geneinfo[rownames(geneinfo) %in% hits,]
tmp2$gene <- rownames(tmp2)
library(ggrepel)
ggplot(tmp2, aes(x = start, xend = stop, y=gene, yend = gene)) +
        geom_segment(linetype = 'solid')+xlab("chrosomose position")+theme_bw()+
        geom_label_repel(aes(label=gene))

# Chr22 ####
ok <- which(geneinfo$chr == "chr22")
aframe <- data.frame(expr[ok,])
means <- do.call(cbind, lapply(infer@observation_grouped_cell_indices, function(x) rowMeans(expr[ok, x])))
means <- data.frame(means)
means$position <- geneinfo$start[ok]
tmp <- reshape2::melt(means, id.vars = "position")
ggplot(tmp, aes(position, value, group = variable, color = variable)) +
        geom_point() +
        geom_hline(yintercept = 1) +
        geom_vline(xintercept = c(18000000,30000000), color = "red") +
        theme_classic()
hits <- rownames(geneinfo)[geneinfo$start < 30000000 & geneinfo$start > 18000000 & geneinfo$chr == "chr22"]
colnames(tmp)[2] <- "sample"
meta <- read.delim("meta/MCL_meta_0118_2022.txt")
meta$sample[meta$sample=="A3_c1"] <- "A3_3prime"
meta$sample[meta$sample=="A3_c2"] <- "A3_5prime"
tmp$clinical_outcome <- meta[match(tmp$sample,meta$sample),'clinical_outcome']
tmp$clinical_outcome <- factor(tmp$clinical_outcome,levels=c("Normal","IBN-S","IBN-Slow","IBN-R","Dual"))
levels(tmp$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
ggplot(tmp[tmp$position < 30000000 & tmp$position > 15000000,], aes(position, value, group = sample, color = clinical_outcome)) +
        geom_point() + ylab("Copy number variation")+
        geom_hline(yintercept = 1) +
        theme_bw()+ggtitle("Chromosome 22")+
        scale_color_manual(values=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))

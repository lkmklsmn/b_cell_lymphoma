setwd("/Users/yanfangfang/Downloads/MW/")
library(Seurat)
library(ggplot2)


# Figure S9B ####
df <- read.delim("trajectory_early_driver_GSEA/Project_ORA_early_driver_2_6_down/enrichment_results_ORA_early_driver_2_6.txt")
df <- df[df$database!="pathway_KEGG",]
df <- df[df$FDR<0.05,c('geneSet','description','enrichmentRatio','FDR')]
df <- df[order(df$enrichmentRatio,decreasing = T),]
colnames(df) <- c('geneSet','description','EnrichmentScore','FDR')
ggplot(df,aes(reorder(description,EnrichmentScore),y=log2(EnrichmentScore))) +
        geom_bar(stat="identity")+ylab("EnrichmentScore")+
        coord_flip()+theme_bw()+xlab("pathway")+
        ggtitle("Dual-R/BTKi-R (trajectory 2/4) vs BTKi-Slow (trajectory 6/7/8), downregulated")+
        theme(axis.title.y =element_blank())

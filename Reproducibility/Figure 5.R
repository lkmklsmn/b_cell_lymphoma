setwd("/Users/yanfangfang/Downloads/MW/")
library(ggplot2)
df1 <- read.delim("inferCNV/infervnv_tumorB_downsample_HMM_outs_include_normal_strigent/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.2.pred_cnv_genes.dat")
df2 <- read.delim("inferCNV/infervnv_tumorB_downsample_HMM_outs_include_normal_strigent/HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.2.pred_cnv_regions.dat")

tmp <- df1[df1$gene %in% c("CCND1","IGH"),]
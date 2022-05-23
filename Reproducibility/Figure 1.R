setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
integrated <- readRDS(file="data/integrated.rds")

# Figure 1A ####
# library("readxl")
# meta <- integrated@meta.data
# meta$sample[meta$sample=="A3" & meta$chemistry=="3prime"] <- "A3_c1"
# meta$sample[meta$sample=="A3" & meta$chemistry=="5prime"] <- "A3_c2"
# df <- meta %>% group_by(chemistry,patient,sample,source,ibrutinib_sensitivity,
#                         days_ibrutinib_treatment,days_ibrutinib_relapse, CART_sensitivity,
#                         days_CART_treatment,celltype) %>%
#         summarise(n = n()) %>% mutate(freq = n / sum(n))
# df <- as.data.frame(df)
# meta$id <- paste0(meta$sample, meta$chemistry, sep = "|")
# df$id <- paste0(df$sample, df$chemistry, sep = "|")
# df$total_no_cell <- table(meta$id)[df$id]
# df$chemistry[df$chemistry=="3prime"] <- "Cohort1"
# df$chemistry[df$chemistry=="5prime"] <- "Cohort2"
# df$clinical_outcome <- factor(df$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
# levels(df$clinical_outcome) <- c("Normal","IBN-S","IBN-Slow","IBN-R","Dual")
# 
# # update days of treatment ####
# meta <- read_excel("meta/Patient info_scRNA_29 Pt_summary_03222021.xlsx",sheet="scRNA_sample info")
# meta <- data.frame(meta[,c(3,5,8,9)])
# colnames(meta) <- c("sample","days_of_ibrutinib_treatment","days_post_ibrutinib_relapse",
#                     "days_post_CAR_T_treatment")
# idx <- match(df$sample,meta$sample)
# df$days_ibrutinib_treatment <- as.numeric(meta[idx,'days_of_ibrutinib_treatment'])
# df$days_post_ibrutinib_relapse <- as.numeric(meta[idx,'days_post_ibrutinib_relapse'])
# df$days_post_CART_treatment <- as.numeric(meta[idx,'days_post_CAR_T_treatment'])
# df$total_no_cell <- as.numeric(df$total_no_cell)
# df$days_ibrutinib_relapse <- df$days_ibrutinib_treatment-df$days_post_ibrutinib_relapse
# df$days_CART_treatment <- df$days_ibrutinib_treatment-df$days_post_CART_treatment
# write.table(df,"meta/MCL_meta_0118_2022.txt",row.names = F,col.names = T,quote=F,sep="\t")
# 
# # create cell type frequency matrix ####
# cell_type_freq <- matrix(NA, length(unique(df$id)), length(unique(df$celltype)))
# rownames(cell_type_freq) <- unique(df$id)
# colnames(cell_type_freq) <- unique(df$celltype)
# tmp <- lapply(unique(df$celltype), function(celltype){
#         lapply(unique(df$id), function(sample){
#                 freq <- unique(df[df$id == sample & df$celltype == celltype, "freq"])
#                 if(length(freq) == 0) freq <- 0
#                 cell_type_freq[sample, celltype] <<- freq
#         })
# })
# 
# # ord <- unique(df$patient[order(df$therapeutic_sensitivity, df$source)])
# # df$patient <- factor(df$patient, levels = ord)
# # 
# # # Calculate max number of days per patient ####
# # asplit <- split(df$days_ibrutinib_treatment, df$patient)
# # max_days <- unlist(lapply(asplit, max))
# # names(max_days) <- names(asplit)
# # df$max_days <- max_days[df$patient]
# 
# # Event ####
# df_uniq <- df[match(unique(df$id), df$id), ]
# event <- do.call(rbind,lapply(unique(df_uniq$patient),function(x){
#         sub <- df_uniq[df_uniq$patient==x,]
#         tmp <- sub[,c("patient","sample","days_ibrutinib_treatment")]
#         if(!(is.na(sub$days_ibrutinib_relapse))){
#                 tmp <- rbind(tmp,c(x,"ibrutinib_relapse",sub$days_ibrutinib_relapse))
#         } else {tmp}
#         if (!(is.na(sub$days_CART_treatment))) {
#                 tmp <- rbind(tmp,c(x,"CART_treatment",sub$days_CART_treatment))
#         } 
#         return(tmp)
# }))
# colnames(event) <- c("patient","event","days")
# event$days <- as.integer(event$days)
# event <- event[!duplicated(event),]
# event[event$event %in% c("A3_c1","A3_c2"),'days'] <- 1005
# 
# # calculate min and max number of days per patient ####
# asplit <- split(event$days, event$patient)
# max_days <- unlist(lapply(asplit, max))
# names(max_days) <- names(asplit)
# event$max_days <- as.integer(max_days[event$patient])
# min_days <- unlist(lapply(asplit, min))
# names(min_days) <- names(asplit)
# event$min_days <- as.integer(min_days[event$patient])
# event$min_days[event$min_days > 0] <- 0
# event$min_days[event$min_days < -50] <- -50
# # add clinical outcome 
# idx <- match(event$event,df$sample)
# event$clinical_outcome <- df[idx,'clinical_outcome']
# 
# idx <- match(event$patient,df$patient)
# event$cohort <- df[idx,'chemistry']
# event$cohort[event$event=="A3_c2"] <- "Cohort2"
# event$group <- "sample"
# event$group[event$event=="ibrutinib_relapse"] <- "ibrutinib_relapse"
# event$group[event$event=="CART_treatment"] <- "CART_treatment"
# event <- event[!(event$event %in% c("Normal_1","Normal_2")),]
# event[event$event=="L6",'days'] <- 981
# event[event$patient=="L",'max_days'] <- 981
# event[event$patient=="Q" & event$event=='ibrutinib_relapse','days'] <- 10
# event[event$patient=="Q" ,'min_days'] <- 0
# event[event$patient=="Q" ,'max_days'] <- 354
# idx <- match(event$patient,df$patient)
# event$cohort <- df[idx,'chemistry']
# event$total_cell <- df[match(event$event,df$sample),'total_no_cell']
# event$total_cell[is.na(event$total_cell)] <- 500
# event$label2 <- event$event
# event$label2[event$group!='sample'] <- NA
# # write.table(event,"meta/event_for_plot_0118_2022.txt",row.names = F,col.names = T,quote=F,sep="\t")

event <- read.delim("meta/event_for_plot_0118_2022.txt")
# if (require(devtools)) {install_github('shabbychef/ggallin')}
event$cohort <- factor(event$cohort)
levels(event$cohort) <- c("Cohort1","Cohort2")
event$cohort[event$event %in% c("L2","L3","L6","A4")] <- "Cohort2"
event$clinical_outcome[event$clinical_outcome=="IBN-S"] <- "BTKi-S"
event$clinical_outcome[event$clinical_outcome=="IBN-Slow"] <- "BTKi-Slow"
event$clinical_outcome[event$clinical_outcome=="IBN-R"] <- "BTKi-R"
event$clinical_outcome[event$clinical_outcome=="Dual"] <- "Dual-R"
event$clinical_outcome <- factor(event$clinical_outcome,levels=c("BTKi-S","BTKi-Slow","BTKi-R","Dual-R"))
library(ggallin)
library(ggrepel)
ggplot(event, aes(days, patient,color=clinical_outcome,label = label2)) +
        scale_x_continuous(trans=pseudolog10_trans,
                           breaks=c(-100,-10,0,10,100,1000))+
        geom_point(aes(shape=group,size=log10(total_cell)))+
        geom_vline(xintercept = 0, linetype = 'dashed') +
        geom_segment(aes(x = min_days, xend = max_days, y = patient, yend = patient), 
                     linetype = 'dashed', color = 'black') + 
        theme_classic()+xlab("Days post BTKi treatment")+
        geom_label_repel(size=3,max.overlaps = 8)+
        facet_wrap(~cohort,ncol = 2,scales = "free")+
        scale_color_manual(values=c("steelblue","lightblue2","peachpuff4","chocolate2"))

# Figure 1C,E ####
integrated$clinical_outcome <- factor(integrated$ibrutinib_sensitivity,levels=c("Normal","S","Slow_responder","R","Dual"))
levels(integrated$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
integrated$clinical_outcome <- factor(integrated$clinical_outcome,
                                      levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R"))
integrated$celltype[integrated$celltype=="Tumor B"] <- "B cells"
DimPlot(integrated,group.by = 'celltype',cols=c("peachpuff4","chocolate3","yellowgreen","wheat3","cyan3",
                                                "khaki","pink"))+
        labs(title="Cell type",subtitle="n=78,740 cells")
DimPlot(integrated,group.by = 'clinical_outcome',
        cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Clinical outcome")

# Figure 1D ####
#features <- c("CCND1","CD79A","CD14","LYZ","CD8A","CD8B","IL7R","LTB","GNLY","HBB","PPBP","FCGR3A","IL3RA","LILRA4")
integrated$celltype <- as.character(integrated$celltype)
integrated$celltype[integrated$celltype=="B cells" & integrated$clinical_outcome!="Normal"] <- "Tumor B"
integrated$celltype[integrated$celltype=="B cells" & integrated$clinical_outcome=="Normal"] <- "Normal B"
Idents(integrated) <- integrated$celltype
m <- FindAllMarkers(integrated,only.pos = T,min.pct = 0.2,min.diff.pct=0.2,
                    logfc.threshold = 0.25,max.cells.per.ident = 200)
top10 <- m %>% group_by(cluster) %>% top_n(3, avg_log2FC)
lst <- top10$gene 
# x1 <- FindMarkers(integrated,only.pos = T,min.pct = 0.2,ident.1 = "Tumor B",
#             ident.2 = "Normal B",logfc.threshold = 0.25,max.cells.per.ident = 200)
lst[1:3] <- c("CCND1","CD79A","CD79B")
# lst[19:21] <- c("CD79A","CD79B","MS4A1")
DotPlot(integrated, features = lst) + RotatedAxis()+
        theme(axis.title.y = element_blank(),axis.title.x=element_blank(),
              axis.text.x = element_text(size = 8,angle = 30, vjust = 0.8, hjust=0.5),
              axis.text.y = element_text(size = 10))+
        scale_color_viridis_c()

# # Figure 1E ####
# library(dplyr)
# meta <- integrated@meta.data
# df <- meta %>% group_by(chemistry,patient,sample,celltype) %>%
#         summarise(n = n()) %>% mutate(freq = n / sum(n))
# df$chemistry[df$chemistry=="3prime"] <- "Cohort1"
# df$chemistry[df$chemistry=="5prime"] <- "Cohort2"
# df <- as.data.frame(df)
# df$celltype <- factor(df$celltype)
# df$celltype <- factor(df$celltype, levels=rev(levels(df$celltype)))
# idx <- match(df$sample,meta$sample)
# df$ibrutinib_sensitivity <- meta[idx,'ibrutinib_sensitivity']
# df$source <- meta[idx,'source']
# df <- arrange(df,factor(source,levels=c("Spleen","LN","BM","PB")))
# lvls <- unique(df$sample)
# df$celltype <- factor(df$celltype,levels=c("pDCs","NK","CD8 T","CD4 T","CD16 Mono","CD14 Mono","B cells"))
# ggplot(df, aes(fill=celltype, y=n, x=factor(sample,levels=lvls))) + 
#         geom_bar(position="stack", stat="identity")+
#         scale_fill_manual(values=c("cyan3", "lightblue4","lightblue3","wheat3",
#                                    "lightblue2","darkgrey","pink"))+
#         ggtitle("Cell type frequencies")+
#         xlab("Sample")+ylab("Number of cells")+theme_bw()+ coord_flip()+
#         facet_wrap(~chemistry,scales = "free")#+
#         #geom_text(aes(label = source),position=position_stack(vjust=2))

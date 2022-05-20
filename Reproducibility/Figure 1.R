# Load R libs ####
library(dplyr)
library(ggallin)
library(ggplot2)
library(ggrepel)
library(Seurat)

# Figure 1A ####
event <- read.delim("../data/event_for_plot_0118_2022.txt")

event$cohort <- factor(event$cohort)
levels(event$cohort) <- c("Cohort1","Cohort2")
event$cohort[event$event %in% c("L2","L3","L6","A4")] <- "Cohort2"
event$clinical_outcome[event$clinical_outcome=="IBN-S"] <- "BTKi-S"
event$clinical_outcome[event$clinical_outcome=="IBN-Slow"] <- "BTKi-Slow"
event$clinical_outcome[event$clinical_outcome=="IBN-R"] <- "BTKi-R"
event$clinical_outcome[event$clinical_outcome=="Dual"] <- "Dual-R"
event$clinical_outcome <- factor(event$clinical_outcome,
                                 levels=c("BTKi-S","BTKi-Slow","BTKi-R","Dual-R"))
ggplot(event, aes(days, patient, color=clinical_outcome, label = label2)) +
        geom_segment(aes(x = min_days, xend = max_days,
                         y = patient, yend = patient), 
                     linetype = 'dashed', color = 'black') + 
        scale_x_continuous(trans=pseudolog10_trans,
                           breaks = c(-100, -10, 0, 10, 100, 1000))+
        geom_point(aes(shape=group, size=log10(total_cell)))+
        geom_vline(xintercept = 0, linetype = 'dashed') +
        theme_classic() + 
        labs(x = "Days post BTKi treatment",
             y = "Patient",
             title = "Longitudinal MCL Patient Sampling") +
        geom_label_repel(size=3, max.overlaps = 10) +
        facet_wrap(~cohort, ncol = 2, scales = "free") +
        scale_color_manual(values=c("steelblue","lightblue2",
                                    "peachpuff4","chocolate2"))



# Figure 1C,E ####
integrated <- readRDS(file="../data/integrated.rds")

integrated$clinical_outcome <- factor(integrated$ibrutinib_sensitivity,
                                      levels=c("Normal","S","Slow_responder","R","Dual"))
levels(integrated$clinical_outcome) <- c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R")
integrated$clinical_outcome <- factor(integrated$clinical_outcome,
                                      levels=c("Normal","BTKi-Fast","BTKi-Slow","BTKi-R","Dual-R"))
integrated$celltype[integrated$celltype=="Tumor B"] <- "B cells"

DimPlot(integrated, 
        group.by = 'celltype',
        cols=c("peachpuff4","chocolate3","yellowgreen","wheat3","cyan3","khaki","pink"))+
        labs(title="Cell type",
             subtitle="n=78,740 cells") +
        theme_void()

DimPlot(integrated,group.by = 'clinical_outcome',
        cols=c("darkgrey","steelblue","lightblue2","peachpuff4","chocolate2"))+
        labs(title="Clinical outcome",
             subtitle="n=78,740 cells") +
        theme_void()

# Figure 1D ####
integrated$celltype <- as.character(integrated$celltype)
integrated$celltype[integrated$celltype=="B cells" &
                            integrated$clinical_outcome!="Normal"] <- "Tumor B"
integrated$celltype[integrated$celltype=="B cells" &
                            integrated$clinical_outcome=="Normal"] <- "Normal B"
Idents(integrated) <- integrated$celltype
m <- FindAllMarkers(integrated,only.pos = T,min.pct = 0.2,min.diff.pct=0.2,
                    logfc.threshold = 0.25,max.cells.per.ident = 200)
top10 <- m %>% group_by(cluster) %>% top_n(3, avg_log2FC)
lst <- top10$gene 
lst[1:3] <- c("CCND1","CD79A","CD79B")

DotPlot(integrated, features = lst) + RotatedAxis()+
        theme(axis.title.y = element_blank(),
              axis.title.x=element_blank(),
              axis.text.x = element_text(size = 10, angle = 90, hjust=1),
              axis.text.y = element_text(size = 10))+
        scale_color_viridis_c()

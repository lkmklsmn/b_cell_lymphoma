setwd("/Users/yanfangfang/Downloads/MW/")
library("Seurat")
library(ggplot2)
library(dplyr)
integrated <- readRDS(file="data/integrated.rds")
integrated <- RunUMAP(integrated,dims=1:10)
integrated <- integrated[,integrated$celltype!="Eryth"]
saveRDS(integrated,"data/integrated.rds")
#mcl <- readRDS(file="data/MCL.rds")
#integrated$celltype_new <- mcl@meta.data[colnames(integrated),'celltype']

# Figure 1A ####
library("readxl")
meta <- integrated@meta.data
df <- meta %>% group_by(chemistry,patient,sample,source,ibrutinib_sensitivity,
                        days_ibrutinib_treatment,days_ibrutinib_relapse, CART_sensitivity,
                        days_CART_treatment,celltype) %>%
        summarise(n = n()) %>% mutate(freq = n / sum(n))
df <- as.data.frame(df)
meta$df <- paste0(meta$sample, meta$chemistry)
df$df <- paste0(df$sample, df$chemistry)
df$total_no_cell <- table(meta$df)[df$df]
df$chemistry[df$chemistry=="3prime"] <- "Cohort1"
df$chemistry[df$chemistry=="5prime"] <- "Cohort2"
df$clinical_outcome <- df$ibrutinib_sensitivity
df$ibrutinib_sensitivity[df$ibrutinib_sensitivity=="Dual"] <- "R"
df$clinical_outcome <- factor(df$clinical_outcome)
levels(df$clinical_outcome) <- c("Dual","Normal","IBN-R","IBN-S","IBN-Slow")
df <- df[,-13]

# update days of treatment ####
meta <- read_excel("Patient info_scRNA_29 Pt_summary_03222021.xlsx",sheet="scRNA_sample info")
meta <- data.frame(meta[,c(3,5,8,9)])
colnames(meta) <- c("sample","days_of_ibrutinib_treatment","days_post_ibrutinib_relapse",
                    "days_of_CAR_T_treatment")
meta$days_of_ibrutinib_treatment <- as.numeric(meta$days_of_ibrutinib_treatment)
idx <- match(df$sample,meta$sample)
df$days_ibrutinib_treatment <- meta[idx,'days_of_ibrutinib_treatment']
df$days_ibrutinib_relapse <- meta[idx,'days_post_ibrutinib_relapse']
df$days_CART_treatment <- meta[idx,'days_of_CAR_T_treatment']

meta <- df
meta$id <- paste(meta$chemistry, meta$sample, sep = "|")
# create cell type frequency matrix ####
cell_type_freq <- matrix(NA, length(unique(meta$id)), length(unique(meta$celltype)))
rownames(cell_type_freq) <- unique(meta$id)
colnames(cell_type_freq) <- unique(meta$celltype)

lapply(unique(meta$celltype), function(celltype){
        lapply(unique(meta$id), function(sample){
                print(paste(celltype, sample))
                freq <- unique(meta[meta$id == sample & meta$celltype == celltype, "freq"])
                if(length(freq) == 0) freq <- NA
                cell_type_freq[sample, celltype] <<- freq
        })
})

# Remove samples with NA values ####
df <- meta[match(unique(meta$id), meta$id), ]
df <- df[!is.na(df$days_ibrutinib_treatment), ]

# ord <- unique(df$patient[order(df$therapeutic_sensitivity, df$source)])
# df$patient <- factor(df$patient, levels = ord)
# 
# # Calculate max number of days per patient ####
# asplit <- split(df$days_ibrutinib_treatment, df$patient)
# max_days <- unlist(lapply(asplit, max))
# names(max_days) <- names(asplit)
# df$max_days <- max_days[df$patient]

# Create plot ####
df$total_no_cell <- as.numeric(df$total_no_cell)
df$days_ibrutinib_relapse <- df$days_ibrutinib_treatment-df$days_ibrutinib_relapse
df$days_CART_treatment <- df$days_ibrutinib_treatment-df$days_CART_treatment
write.table(df,"MCL_meta.txt",row.names = F,col.names = T,quote=F,sep="\t")

event <- do.call(rbind,lapply(unique(df$patient),function(x){
        sub <- df[df$patient==x,]
        tmp <- sub[,c("patient","sample","days_ibrutinib_treatment")]
        if(!(is.na(sub$days_ibrutinib_relapse))){
                tmp <- rbind(tmp,c(x,"ibrutinib_relapse",sub$days_ibrutinib_relapse))
        } else {tmp}
        if (!(is.na(sub$days_CART_treatment))) {
                tmp <- rbind(tmp,c(x,"CART_treatment",sub$days_CART_treatment))
        } 
        return(tmp)
}))
colnames(event) <- c("patient","event","days")
write.table(event,"event.txt",row.names = F,col.names = T,quote=F,sep="\t")

library(ggrepel)
ggplot(df, aes(days_ibrutinib_treatment, patient,color = clinical_outcome, label = sample)) +
        geom_point(aes(size=total_no_cell))+
        facet_wrap(~ chemistry, ncol = 2, scales = "free_y") +
        geom_vline(xintercept = 0, linetype = 'dashed') +
        #geom_segment(aes(x = -203, xend = max_days, y = patient, yend = patient, size = 1), linetype = 'dashed', color = 'black') + 
        theme_classic()+xlab("Days post Ibrutinib treatment")+
        geom_label_repel()
# Create pie chart ####
data <- data.frame(group=LETTERS[1:5],value=c(13,7,9,21,2))

# Basic piechart
data <- data.frame(value = cell_type_freq[1,],group = colnames(cell_type_freq))
data$value[is.na(data$value)] <- 0

ggplot(data, aes(x="", y=value, fill=group)) +
        geom_bar(stat="identity", width=1) +
        coord_polar("y", start=0)

library(scatterpie)
long <- rnorm(50, sd=100)
lat <- rnorm(50, sd=50)
d <- data.frame(long=long, lat=lat)
d <- with(d, d[abs(long) < 150 & abs(lat) < 70,])
n <- nrow(d)
d$region <- factor(1:n)
d$A <- abs(rnorm(n, sd=1))
d$B <- abs(rnorm(n, sd=2))
d$C <- abs(rnorm(n, sd=3))
d$D <- abs(rnorm(n, sd=4))
d[1, 4:7] <- d[1, 4:7] * 3
head(d)

cell_type_freq[is.na(cell_type_freq)] <- 0
final <- data.frame(df, cell_type_freq[df$id,])

final$patient_num <- match(final$patient, unique(final$patient))

ggplot(data= final, aes(x=days_ibrutinib_treatment, patient_num)) +
        geom_scatterpie(data = final, aes(x=days_ibrutinib_treatment, patient_num, r = log10(total_no_cell)),
                        cols=gsub(" ", ".", fixed = T, colnames(cell_type_freq))) + coord_


# Figure 1B,C,D ####
integrated$ibrutinib_sensitivity <- factor(integrated$ibrutinib_sensitivity)
levels(integrated$ibrutinib_sensitivity) <- c("Dual","Normal","IBN-R","IBN-S","IBN-Slow")
DimPlot(integrated,group.by = 'celltype',cols=c("chocolate2","yellowgreen","wheat3","cyan3",
                                                "khaki","pink","peachpuff4"))
DimPlot(integrated,group.by = 'ibrutinib_sensitivity',
        cols=c("coral4","grey","tomato3","peachpuff4","skyblue4"))+
        labs(title="Clinical outcomes")
features <- c("CCND1","CD79A","CD14","LYZ","CD8A","CD8B","IL7R","LTB","GNLY","HBB","PPBP","FCGR3A","IL3RA","LILRA4")
DotPlot(integrated, features = features) + RotatedAxis()

# Figure 1E ####
library(dplyr)
meta <- integrated@meta.data
df <- meta %>% group_by(chemistry,patient,sample,celltype) %>%
        summarise(n = n()) %>% mutate(freq = n / sum(n))
df$chemistry[df$chemistry=="3prime"] <- "Cohort1"
df$chemistry[df$chemistry=="5prime"] <- "Cohort2"
df <- as.data.frame(df)
df$celltype <- factor(df$celltype)
df$celltype <- factor(df$celltype, levels=rev(levels(df$celltype)))
idx <- match(df$sample,meta$sample)
df$ibrutinib_sensitivity <- meta[idx,'ibrutinib_sensitivity']
df$source <- meta[idx,'source']
df <- arrange(df,factor(source,levels=c("Spleen","LN","BM","PB")))
lvls <- unique(df$sample)
df$celltype <- factor(df$celltype,levels=c("pDCs","NK","Eryth","CD8 T","CD4 T","CD16 Mono","CD14 Mono","Tumor B"))
ggplot(df, aes(fill=celltype, y=n, x=factor(sample,levels=lvls))) + 
        geom_bar(position="stack", stat="identity")+
        scale_fill_manual(values=c("cyan3", "lightblue4","lightblue3","wheat3",
                                   "lightblue2","darkgrey","peachpuff","pink"))+
        ggtitle("Cell type proportations")+
        xlab("Sample")+ylab("Number of cells")+theme_bw()+ coord_flip()+
        facet_wrap(~chemistry,scales = "free")

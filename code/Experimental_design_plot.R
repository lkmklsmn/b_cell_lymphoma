# Load R libs ####
library(readxl)
library(ggplot2)

# Load metadata ####
setwd("/Users/lukas/OneDrive/Documents/GitHub/b_cell_lymphoma/")
meta <- read.csv("data/meta_new.csv")

meta$id <- paste(meta$chemistry, meta$sample, sep = "|")
meta$days_ibrutinib_treatment <- as.numeric(meta$days_ibrutinib_treatment)

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
tmp <- meta[match(unique(meta$id), meta$id), ]
tmp <- tmp[!is.na(tmp$days_ibrutinib_treatment), ]

# Sort patients ####
ord <- unique(tmp$patient[order(tmp$therapeutic_sensitivity, tmp$source)])
tmp$patient <- factor(tmp$patient, levels = ord)

# Calculate max number of days per patient ####
asplit <- split(tmp$days_ibrutinib_treatment, tmp$patient)
max_days <- unlist(lapply(asplit, max))
names(max_days) <- names(asplit)

tmp$max_days <- max_days[tmp$patient]

# Create plot ####
p <- ggplot(tmp, aes(days_ibrutinib_treatment, patient,
                 color = therapeutic_sensitivity, label = sample, size = log10(total_no_cell))) +
  facet_wrap(~ chemistry, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_segment(aes(x = -203, xend = max_days, y = patient, yend = patient, size = 1), linetype = 'dashed', color = 'black') + 
  geom_label() +
  theme_classic()

# Create pie chart ####
data <- data.frame(
  group=LETTERS[1:5],
  value=c(13,7,9,21,2)
)

# Basic piechart
data <- data.frame(value = cell_type_freq[1,],
                   group = colnames(cell_type_freq))
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
final <- data.frame(tmp, cell_type_freq[tmp$id,])

final$patient_num <- match(final$patient, unique(final$patient))

ggplot(data= final, aes(x=days_ibrutinib_treatment, patient_num)) +
  geom_scatterpie(data = final, aes(x=days_ibrutinib_treatment, patient_num, r = log10(total_no_cell)),
                  cols=gsub(" ", ".", fixed = T, colnames(cell_type_freq))) + coord_



setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
options(scipen= 99)
require(pheatmap)
require(data.table)

#----------------------------------------------#
# 1- Correlations counts
#----------------------------------------------#
counts <- data.table(file= list.files("db/counts", ".txt", recursive = T, full.names = T))
counts <- counts[grepl("epiCancer", file)]
counts[, cdition:= gsub(".txt", "", basename(file))]
counts <- counts[, fread(file, col.names = c("rn", "counts")), c(colnames(counts))]
counts <- counts[!grepl("^__", rn)]
dmat <- dcast(counts, rn~cdition, value.var = "counts")
mat <- as.matrix(dmat, 1)
cor <- cor(mat)

pheatmap(cor, filename =  "pdf/correlations_transcriptomes_replicates_heatmap.pdf", width = 9, height = 8, display_numbers = T)






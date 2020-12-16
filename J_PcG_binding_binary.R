setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)
require(pheatmap)
require(rtracklayer)
require(readxl)

peaks <- as.data.table(read_xlsx("../../public_data/dm6/peaks/PcG_peaks_Loubiere_SA_2020.xlsx", skip = 2))

# Import
dat <- readRDS("Rdata/som_clustering_transcriptomes.rds")

t1 <- fread("db/bam/PH29_rep3.bam.junction.bed", header= F)
t2 <- fread("db/bam/W29_rep1.bam.junction.bed", header= F)

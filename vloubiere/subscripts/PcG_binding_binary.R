peaks <- as.data.table(read_xlsx("../../public_data/dm6/peaks/PcG_peaks_Loubiere_SA_2020.xlsx", skip = 2))

# Import
dat <- readRDS("Rdata/som_clustering_transcriptomes.rds")

t1 <- fread("db/bam/PH29_rep3.bam.junction.bed", header= F)
t2 <- fread("db/bam/W29_rep1.bam.junction.bed", header= F)

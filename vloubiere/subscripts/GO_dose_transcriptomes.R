setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")


dat <- readRDS("Rdata/clustering_dose_transcriptomes.rds")
.l <- split(dat$FBgn, dat$rcl)
.l <- lapply(.l, unique)

pdf("pdf/clustering/clustering_dose_GOs.pdf", height = 10, width = 10)
vl_GO_clusters(.l, 
               cex = 0.3, 
               N_top = 10)
dev.off()
setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(cdition, system, FC_file)]
dat[, cdition:= factor(cdition, 
                       levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dat <- dat[diff!="unaffected"]

pdf("pdf/Figures/upsetplot_overlapping_genes_RNA.pdf", width = 7.5)
dat[, {
  vl_upset_plot(split(FBgn, list(cdition, diff)), 
                intersection_cutoff = 10)
}, system]
dev.off()

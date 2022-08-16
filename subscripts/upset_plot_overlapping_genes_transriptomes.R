setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA"]
meta <- na.omit(meta[, .(cdition= gsub("^RNA_", "", cdition), FC_file)])
meta[, cdition:= factor(cdition, 
                        levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dat <- meta[, fread(FC_file), (meta)]
dat <- dat[diff!="unaffected"]

pdf("pdf/Figures/upsetplot_overlapping_genes_RNA.pdf", width = 7.5)
vl_upset_plot(split(dat$FBgn, dat[, .(cdition, diff)]), 
              intersection_cutoff = 10)
dev.off()

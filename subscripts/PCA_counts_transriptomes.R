setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import data
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[, as.data.table(counts(readRDS(dds_file))), .(dds_file, system)]
dat <- dat[, {
  mat <- t(log2(.SD+1))
  as.data.table(prcomp(mat)$x, keep.rownames = "cdition")
}, .(dds_file, system)]
dat[, cdition:= gsub("_rep1|_rep2|_rep3", "", cdition)]
dat[, col:= rainbow(.NGRP)[.GRP], cdition]

pdf("pdf/Figures/PCA_counts_RNA.pdf",
    width = 5,
    height = 5.5)
par(las= 1)
dat[, {
  plot(PC1, 
       PC2, 
       col= adjustcolor(col, 0.6), 
       pch= 16,
       cex= 2,
       main= system)
  legend(c("bottomright", "topright")[.GRP],
         unique(cdition),
         col= adjustcolor(unique(col), 0.6),
         pch= 16,
         bty= "n",
         pt.cex= 2)
  print("")
}, system]
dev.off()
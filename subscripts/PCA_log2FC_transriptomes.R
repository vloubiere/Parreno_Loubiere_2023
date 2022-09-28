setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import data
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(cdition, FC_file, dds_file, system)]
dat[, cdition:= paste0(cdition, "_", system)]
mat <- as.matrix(dcast(dat, FBgn~cdition, value.var = "log2FoldChange"), 1)
mat <- na.omit(mat)
mat <- as.data.table(prcomp(t(mat))$x, keep.rownames = "cdition")
mat[, c("cdition", "system"):= tstrsplit(cdition, "_")]
mat[, col:= vl_palette_few_categ(.NGRP)[.GRP], cdition]
mat[, pch:= c(16, 15)[.GRP], system]

pdf("pdf/RNA_PCA_log2FC.pdf",
    width = 5,
    height = 5.5)
par(las= 1)
mat[, {
  plot(PC1, 
       PC2, 
       col= adjustcolor(col, 0.6), 
       pch= pch,
       cex= 2,
       main= "PCA compare GFP/noGFP systems")
  legend("topleft",
         paste0(cdition, "_", system),
         col= adjustcolor(col, 0.6),
         pch= pch,
         bty= "n",
         pt.cex= 1,
         cex= 0.8)
  print("")
}]
dev.off()

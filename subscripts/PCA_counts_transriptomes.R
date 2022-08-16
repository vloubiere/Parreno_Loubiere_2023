setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import data
meta <- fread("Rdata/processed_metadata_RNA.txt")
dds <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA", unique(dds_file)]
mat <- counts(readRDS(dds))
mat <- t(log2(mat+1))
pca <- as.data.table(prcomp(mat)$x, keep.rownames = "cdition")
pca[, cdition:= tstrsplit(cdition, "_", keep= 1)]
pca[, col:= vl_palette_many_categ(.NGRP)[.GRP], cdition]

pdf("pdf/Figures/PCA_counts_RNA.pdf",
    width = 5,
    height = 5.5)
par(las= 1)
pca[, {
  plot(PC1, 
       PC2, 
       col= adjustcolor(col, 0.6), 
       pch= 16,
       cex= 2)
  legend("topright",
         unique(cdition),
         col= adjustcolor(unique(col), 0.6),
         pch= 16,
         bty= "n",
         pt.cex= 2)
}]
dev.off()
setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
pdf("pdf/Figures/PCA_FPKMs_RNA.pdf",
    width = 5,
    height = 5.5)
for(geno in c("GFP-", "GFP+"))
{
  dds <- switch(geno,
                "GFP-"= "db/dds/RNA/epiCancer_ED_GFP-_system_RNA_dds.rds",
                "GFP+"= "db/dds/RNA/epiCancer_ED_GFP+_system_RNA_dds.rds")
  dat <- DESeq2::fpkm(readRDS(dds))
  dat <- melt(as.data.table(dat, keep.rownames = T), id.vars = "rn")
  dat[, variable:= gsub("_1.bam|_2.bam|_3.bam", "", variable)]
  dat <- dat[, .(FPKM= mean(value)), .(rn, variable)]
  dat[, variable:= gsub("^TPH", "PH", variable)]
  dat[, variable:= gsub("^TW", "W", variable)]
  dat[variable=="PHJ9", variable:= "PHD9"]
  dat[variable=="PHJ11", variable:= "PHD11"]
  dat[, col:= vl_palette_many_categ(.NGRP)[.GRP], variable]
  Cc <- unique(dat[, .(variable, col)])
  dat <- dcast(dat, rn~variable, value.var = "FPKM")
  mat <- as.matrix(dat, 1)
  mat[is.na(mat)] <- 0
  mat <- scale(mat)
  pca <- prcomp(mat)$rotation
  
  #Plot
  plot(pca[,"PC1"],
       pca[,"PC2"],
       col= adjustcolor(Cc[rownames(pca), col, on= "variable"], 0.7),
       pch= ifelse(grepl("GFP\\+", rownames(pca)), 15, 19),
       las= 1,
       cex= 1.5,
       main= geno)
  legend("topleft", 
         rownames(pca),
         pch= ifelse(grepl("GFP\\+", rownames(pca)), 15, 19),
         col= adjustcolor(Cc[rownames(pca), col, on= "variable"], 0.7),
         bty= "n",
         cex= 0.7)
}
dev.off()
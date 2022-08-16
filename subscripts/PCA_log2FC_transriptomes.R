setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object %in% c("epiCancer_ED_GFP-_system_RNA", "RNA_phRNAi_SA2020_ED_Loubiere")]
meta <- na.omit(meta[, .(cdition= gsub("^RNA_", "", cdition), FC_file)])
meta <- unique(meta)
meta[cdition=="PHRNAI_ED", cdition:= "PH25"]
dat <- meta[, fread(FC_file), (meta)]
mat <- dcast(dat, 
             FBgn~cdition, 
             value.var = "log2FoldChange")
mat <- t(as.matrix(na.omit(mat), 1))
pca <- as.data.table(prcomp(mat)$x, keep.rownames = "cdition")
pca[, col:= vl_palette_few_categ(.NGRP)[.GRP], cdition]
pca <- pca[order(factor(cdition, c("PH18", "PHD9", "PHD11", "PH29", "PH25")))]


pdf("pdf/Figures/PCA_log2FC_RNA.pdf",
    width = 5,
    height = 5.5)
par(las= 1)
pca[, {
  plot(PC1, 
       PC2, 
       col= adjustcolor(col, 0.6), 
       pch= 16,
       cex= 2,
       ylim= c(-40, 50))
  legend("topleft",
         unique(cdition),
         col= adjustcolor(unique(col), 0.6),
         pch= 16,
         bty= "n",
         pt.cex= 2)
}]
pca[cdition!="PH25", {
  plot(PC1, 
       PC2, 
       col= adjustcolor(col, 0.6), 
       pch= 16,
       cex= 2,
       ylim= c(-40, 45))
  legend("topleft",
         unique(cdition),
         col= adjustcolor(unique(col), 0.6),
         pch= 16,
         bty= "n",
         pt.cex= 2)
}]
dev.off()

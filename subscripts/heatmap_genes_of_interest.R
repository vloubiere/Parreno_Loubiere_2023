setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- readRDS("Rdata/final_gene_features_table.rds")

# Import genes of interest ----
genes <- fread("Rdata/list_genes_interest_heatmaps.txt", fill= T)
genes[, cluster:= group]
genes[group %in% c("PRC1", "PRC2"), group:= "PcG"]
genes[, group:= factor(group, c("PcG", "HOX", "Eye Disc Dev.", "JNK pathway"))]

# Merge ----
dat <- merge(genes,
             dat[, .(FBgn, log2FoldChange_PH18, log2FoldChange_PH29, log2FoldChange_PHD9, log2FoldChange_PHD11)],
             all.x= T)
setnames(dat,
         c("log2FoldChange_PH18", "log2FoldChange_PH29", "log2FoldChange_PHD9", "log2FoldChange_PHD11"),
         c("No ph-KD", "Constant ph-KD", "Transient ph-KD (d9)", "Transient ph-KD (d11)"))

# Order ----
setorderv(dat,
          c("group", "symbol"))

# Plot ----
pdf("pdf/review_heatmaps_genes_of_interest.pdf",
    height= 3,
    width = 5.8)
vl_par(bty= "o",
       mai= c(.9,.4,.2,.1),
       omi= c(0,0,0,1),
       mfrow= c(1,4),
       font.main= 1)
dat[, {
  mat <- as.matrix(.SD, 1)
  vl_heatmap(mat,
             tilt.colnames= T,
             breaks= c(-4,0,4),
             display.numbers= T,
             display.numbers.matrix= round(mat, 1),
             display.numbers.cex= .5,
             cluster.cols= F,
             grid= F,
             box.lwd = 0.75,
             row.clusters= if(group=="PcG") cluster else 1,
             show.legend= .GRP==.NGRP,
             legend.title= "FoldChange (log2)",
             legend.cex= .7,
             main= as.character(group))
  .SD
}, group,
.SDcols= c("symbol", "No ph-KD", "Constant ph-KD", "Transient ph-KD (d9)", "Transient ph-KD (d11)")]
dev.off()

setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data ----
dat <- readRDS("Rdata/final_ATAC_table.rds")
dat <- dat[cl!="Unaffected"]
setorderv(dat, "cl")
dat[, cl:= paste0(cl, " (", formatC(.N, big.mark = ","), ")"), cl]
dat[, cl:= droplevels(cl)]

# Plot ----
pdf("pdf/review_ATAC_clustering.pdf", 4, 3)
vl_par(mai= c(1.2,1.7,.6,1.7),
       las= 2,
       mgp= c(0,0,0),
       cex.axis= 9/12)
vl_heatmap(dat[, .("Constant ph-KD"= log2FoldChange_PH29,
                   "Transient ph-KD"= log2FoldChange_PHD11)], 
           row.clusters = dat$cl,
           cluster.rows= F,
           cluster.cols= F, 
           breaks = c(-2.5, -0.25, 0.25, 2.5), 
           col= c("cornflowerblue", "white", "white", "tomato"),
           show.rownames = F, 
           legend.title = "Fold Change (log2)",
           legend.cex = 6/12,
           show.col.clusters = F,
           tilt.colnames = T,
           row.clusters.pos = "left")
dev.off()
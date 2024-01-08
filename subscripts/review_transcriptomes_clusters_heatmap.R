setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import clustering data and extra features ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[cl!="Unaffected"]
setorderv(dat, "cl")
dat[, cl:= paste0(cl, " (", .N, ")"), cl]
dat[, cl:= droplevels(cl)]

# PLOT ----
pdf("pdf/Figure_2_clustering.pdf", 4.6, 3.5)
# Heatmap
par(mai= c(1.3,2,.2,1.7),
    las= 2,
    mgp= c(0,0,0))
vl_heatmap(dat[, .(log2FoldChange_PH29, log2FoldChange_PHD9, log2FoldChange_PHD11)],
           row.clusters = dat$cl,
           cluster.rows= F,
           cluster.cols= F,
           breaks = c(-4, 0, 4),
           show.rownames = F,
           tilt.colnames= T,
           col = c("cornflowerblue", "white", "tomato"),
           legend.title = "Fold Change (log2)",
           legend.cex = .7, 
           row.clusters.pos = "left")
dev.off()
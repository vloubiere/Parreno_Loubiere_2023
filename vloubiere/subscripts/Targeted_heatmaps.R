dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")

Cc <- c("blue", "cornflowerblue", "white", "tomato", "red")
lim <- c(-5,0,5)

dat[, pdf:= tempfile(fileext = "pdf"), cl]
dat[, {
  n_genes <- length(unique(FBgn))
  
  pdf(pdf, height = 0.275*n_genes+2)
  par(mar=c(8,8,2,2))
  mat <- as.matrix(dcast(.SD, symbol~cdition, value.var = "log2FoldChange"), 1)
  vl_heatmap(mat, 
             cluster_rows = F,
             cluster_cols = F, 
             col= Cc, 
             breaks = lim,
             legend_title = "log2FC", 
             display_numbers = T)
  mtext(paste0("Cluster ", cl, " (", n_genes, " genes)"))
  dev.off()
}, .(cl, pdf)]

pdf_combine(unique(dat$pdf), output = "pdf/clusters_targeted_heatmaps.pdf")

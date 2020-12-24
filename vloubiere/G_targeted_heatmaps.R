setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)
require(pheatmap)
require(ontologyIndex)
require(STRINGdb)
require(pdftools)

dat <- readRDS("Rdata/som_clustering_transcriptomes.rds")

Cc <- c("cornflowerblue", "white", "tomato")
lim <- c(-5,0,5)

dat[, pdf:= tempfile(fileext = "pdf"), cl]
dat[, {
  n_genes <- length(unique(FBgn))
  pdf(pdf, height = 0.275*n_genes+2)
  res <- my_heatmap(.SD, row.BY = "symbol", col.BY = "cdition", value.var = "log2FoldChange", col = Cc, cluster_cols = F, 
                    leg_title = "log2FC", breaks = lim, main = paste0("Cluster ", cl, " (", n_genes, " genes)"))
  text(res$xcoor, res$ycoor, round(res$log2FoldChange, 1))
  dev.off()
}, .(cl, pdf)]

pdf_combine(unique(dat$pdf), output = "pdf/clusters_targeted_heatmaps.pdf")

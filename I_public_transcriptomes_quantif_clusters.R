setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)
require(pheatmap)
require(ontologyIndex)
require(rtracklayer)
require(digest)

dat <- readRDS("Rdata/som_clustering_transcriptomes.rds")
dat <- unique(dat[, .(FBgn, cl, symbol, ycoor)])
dat <- dat[, .(file= c("db/FC_tables/RNA_development_RNA_72hED_vs_RNA_WTE1416_FC.txt",
                       "db/FC_tables/RNA_development_RNA_96hED_vs_RNA_72hED_FC.txt",
                       "db/FC_tables/RNA_development_RNA_120hED_vs_RNA_96hED_FC.txt")), dat]
dat <- dat[, { 
  .c <- fread(file)
  cbind(.SD, .c[match(FBgn, V1), baseMean:padj])
}, file]
dat[, cdition := gsub("RNA_development_|_FC.txt|_RNA|RNA_", "", basename(file))]
dat[, cdition := gsub("_vs_", "/", cdition)]

setorderv(dat, "ycoor", -1)
my_heatmap(dat, row.BY = "ycoor", col.BY = "cdition", value.var = "log2FoldChange", cluster_rows = F, cluster_cols = F, breaks = c(-5, 0, 5))
abline(h= dat[, max(ycoor), cl]$V1)

pdf("pdf/dev_transcriptomes_clusters.pdf", width = 5)
mat <- as.matrix(dcast(dat, cl~cdition, value.var = "log2FoldChange", fun.aggregate = mean, na.rm= T), 1)
mat <- mat[, c(2,3,1)]
pl <- my_heatmap(mat, leg_title = "mean zscore", cluster_rows = F, cluster_cols = F, breaks = c(-3,0,3))
text(pl$xcoor, pl$ycoor, round(pl$plot_var, 1))
dev.off()

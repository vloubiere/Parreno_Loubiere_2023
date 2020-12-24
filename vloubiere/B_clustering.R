setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)
require(pheatmap)
require(rtracklayer)

# Import
dat <- data.table(files= sapply(c("EzJ9_vs_WKD", "EzJ11_vs_WKD", "Ez29_vs_W29", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29"), function(x) 
  list.files("db/FC_tables/", x, full.names = T)))
dat[, cdition:= tstrsplit(files, "epiCancer_|_FC", keep= 2)]
dat <- dat[, fread(files, select = c(1,2,3,6), col.names = c("FBgn", "baseMean", "log2FoldChange", "padj")),  (dat)]
# Remove unchanged genes
sel <- dat[, any(padj<0.001 & abs(log2FoldChange)>log2(1.5)), FBgn][(V1), FBgn]
dat <- dat[FBgn %in% sel]
# Remove genes deferentially expressed at 18
sel <- rbind(fread("/_R_data/projects/epigenetic_cancer/db/FC_tables/RNA_epiCancer_PH18_vs_W18_FC.txt"),
             fread("/_R_data/projects/epigenetic_cancer/db/FC_tables/RNA_epiCancer_Ez18_vs_W18_FC.txt"))
sel <- sel[padj>=0.01, V1]
dat <- dat[FBgn %in% sel]
# Clip extremes
dat[, clipped_log2FoldChange:= log2FoldChange]
dat[log2FoldChange> 8, clipped_log2FoldChange:= 8]
dat[log2FoldChange< -8, clipped_log2FoldChange:= -8]
# Clustering
grid <- somgrid(4, 4, "hexagonal", toroidal= T)
train <- as.matrix(dcast(dat, FBgn~cdition, value.var = "clipped_log2FoldChange"), 1)
set.seed(5)
som <- supersom(train, grid= grid, maxNA.fraction = 0.1, user.weights = c(1, 10))
# Add to dat
dat[, cl:= som$unit.classif[match(dat$FBgn, rownames(som$data[[1]]))]]
setorderv(dat, "cl")
# Add genes symbols
symbols <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
dat[symbols, symbol:= i.gene_symbol, on= "FBgn==gene_id "]
# Aggregate
agg <- dat[, .(mean_log2FC= mean(clipped_log2FoldChange), N= .N), .(cdition, cl= paste0("cluster ", cl, " ("))]
agg[, cl:= paste0(cl, N, ")")]

pdf("pdf/clustering_transcriptomes.pdf", width = 10)
par(mfrow= c(1,2), mar= c(8.1,8.1,5.1,6.1))
Cc <- c("blue", "cornflowerblue", "white", "tomato", "red")
res <- my_heatmap(dat, row.BY = "FBgn", col.BY = "cdition", value.var = "clipped_log2FoldChange", col = Cc, cluster_rows = F, row_labels = F, 
                  cluster_cols = F, leg_title = "log2FC")
cl_pos <- rev(dat[, length(unique(FBgn)), cl]$V1)
abline(h= cumsum(cl_pos), lwd= 0.5)
axis(2, at = cumsum(cl_pos)-cl_pos/2, rev(unique(dat$cl)), las= 1, tick = 0, cex.lab= 0.2)
my_heatmap(agg, row.BY = "cl", col.BY = "cdition", value.var = "mean_log2FC", col = Cc[2:4], cluster_rows = F, breaks = c(-5, 0, 5), 
           cluster_cols = F, leg_title = "log2FC")
dev.off()

saveRDS(res, "Rdata/som_clustering_transcriptomes.rds")

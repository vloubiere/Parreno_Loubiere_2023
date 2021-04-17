dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")


Cc <- c("blue", "cornflowerblue", "white", "tomato", "red")
my_heatmap(dat, 
           row.BY = "FBgn", 
           col.BY = "cdition", 
           value.var = "clipped_log2FoldChange", 
           col = Cc, 
           cluster_rows = F, 
           row_labels = F, 
           cluster_cols = F, 
           leg_title = "log2FC")
          
dat <-
# Add PRC1 binding
PRC1 <- fread("Rdata/PcG_binding_assignment_final.txt")
TSS <- PRC1[grepl("^TSS", ID) & !is.na(PRC1_binding)]
clean$PRC1_prom_binding <- log2(TSS[clean, .N, .EACHI, on= "gene_id==FBgn"]$N+1)
RE <- PRC1[grepl("^RE", ID) & !is.na(PRC1_binding)]
clean$PRC1_RE_binding <- log2(RE[clean, .N, .EACHI, on= "gene_id==FBgn"]$N+1)


# Add PRC1 binding
mat <- as.matrix(unique(res[cdition==cdition[1], .(ycoor, PRC1_prom_binding, PRC1_RE_binding)]), 1)
mat[mat>1.5] <- 1.5
my_heatmap(mat, 
           col = c("blue", "yellow"),
           breaks = c(0, 1.5), 
           cluster_rows = F, 
           row_labels = F, 
           cluster_cols = F,
           leg_title = "PRC1 binding")
abline(h= cumsum(cl_pos), lwd= 0.5)
axis(2, at = cumsum(cl_pos)-cl_pos/2, rev(unique(clean$cl)), las= 1, tick = 0, cex.lab= 0.2)


# Add PRC1 binding
mat <- as.matrix(unique(r2[cdition==cdition[1], .(ycoor, perc_prom_bound, perc_RE_bound)]), 1)
my_heatmap(mat, 
           col = c("blue", "yellow"),
           cluster_rows = F, 
           row_labels = F, 
           cluster_cols = F,
           leg_title = "PRC1 binding")

#### ChIP


setorderv(.q, "ycoor", -1)
my_heatmap(.q, row.BY = "ycoor", col.BY = "bw_file", value.var = "norm_counts", cluster_rows = F, breaks = c(-2,0, 2))
abline(h= .q[, max(ycoor), cl]$V1)

pdf("pdf/HTM_clusters.pdf", width = 8)
mat <- as.matrix(dcast(.q, cl~bw_file, value.var = "bw_counts", fun.aggregate = mean), 1)
my_heatmap(scale(mat), cluster_rows = F, leg_title = "mean zscore")
dev.off()


# Transcriptomes
setorderv(dat, "ycoor", -1)
my_heatmap(dat, row.BY = "ycoor", col.BY = "cdition", value.var = "log2FoldChange", cluster_rows = F, cluster_cols = F, breaks = c(-5, 0, 5))
abline(h= dat[, max(ycoor), cl]$V1)

pdf("pdf/dev_transcriptomes_clusters.pdf", width = 5)
mat <- as.matrix(dcast(dat, cl~cdition, value.var = "log2FoldChange", fun.aggregate = mean, na.rm= T), 1)
mat <- mat[, c(2,3,1)]
pl <- my_heatmap(mat, leg_title = "mean zscore", cluster_rows = F, cluster_cols = F, breaks = c(-3,0,3))
text(pl$xcoor, pl$ycoor, round(pl$plot_var, 1))
dev.off()


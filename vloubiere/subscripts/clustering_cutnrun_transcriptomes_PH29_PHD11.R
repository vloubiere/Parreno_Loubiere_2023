setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)

# Import gene symbols
if(!exists("FBGN"))
{
  FBGN <- as.data.table(rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf"))
  FBGN <- unique(FBGN[, .(gene_id, gene_symbol)])
  setkeyv(FBGN, "gene_id")
}

# Import transcriptome data
dat <- readRDS("Rdata/RNA_tables_object.rds")$FC$cutnrun_genotype
# REmove genes affected in PH18
dat <- dat[!FBgn %in% dat[cdition=="PH18" & padj<0.01, FBgn]]
# Select conditions to cluster
dat <- dat[cdition %in% c("PH29", "PHD11")]
# Use only misregulated genes
dat <- dat[, .(check= any(padj<0.05 & abs(log2FoldChange)>1), 
               cdition, 
               log2FoldChange, 
               padj), FBgn][(check)]
# dcast FC
mat <- dcast(dat, FBgn~cdition, value.var = "log2FoldChange")
mat <- as.matrix(mat, 1)
# Clip outliers
mat[mat>10] <- 8
mat[mat<(-10)] <- -8
# dcast padj
pmat <- dcast(dat, FBgn~cdition, value.var = "padj")
pmat <- as.matrix(pmat, 1)
# Use pmat to make bin table
bin <- pmat
bin[] <- 0
bin[pmat<0.05 & apply(mat, 2, between, 1, Inf)] <- 1 #sig up
bin[pmat<0.05 & apply(mat, 2, between, -Inf, -1)] <- -1 #sig down
# Clustering
seed <- 2
grid <- somgrid(5,5,"hexagonal",toroidal= T)
set.seed(seed)
som <- supersom(data= list(bin,
                           mat), 
                grid= grid, 
                user.weights= c(10,2))
Nhc <- 7
hc <- cutree(hclust(dist(som$codes[[1]])), Nhc)[som$unit.classif]
# Interaction network
Cc <- igraph::categorical_pal(Nhc)
net <- vl_STRING_interaction(symbols = FBGN[rownames(mat), gene_symbol],
                             size = abs(rowMeans(mat)),
                             col= Cc[hc],
                             score_cutoff = 900)
# SPlit diff expressed genes
diff <- dat[padj<0.05 & abs(log2FoldChange)>1]
diff <- split(diff$FBgn, paste0(diff$cdition, ifelse(diff$log2FoldChange>0, " UP", " DOWN")))
# SAVE
obj <- cbind(mat, pmat)
colnames(obj) <- paste0(colnames(obj), rep(c("_log2FoldChange", "_padj"), each= 2))
obj <- as.data.table(obj, keep.rownames = "FBgn", key= "FBgn")
obj[, symbol:= FBGN[FBgn, gene_symbol]]
obj[, cl:= hc]
saveRDS(obj,
        "Rdata/clustering_cutnrun_genotype_transcriptomes_PH29_PHD11.rds")

#------------------------------------------------#
# PLOT
#------------------------------------------------#
pdf("pdf/clustering/clustering_cutnrun_genotype_transcriptomes_PHD11_PH29.pdf", width = 16)
layout(matrix(c(1,1,1,2,3,3,3,3,3,3,4,4,4,4,4), nrow= 1))
# upset plot
vl_upset_plot(diff)
# Clustering heatmap
vl_heatmap(mat[order(hc),],
           cluster_rows = F,
           breaks = c(-5, -2, -0.5, 0, 0.5, 2, 5),
           col= c("blue", "cornflowerblue", "white", "white", "white", "tomato", "red"), 
           show_rownames = F)
abline(h= nrow(mat)-cumsum(table(hc))+0.5)
# GO
par(cex= 0.7)
vl_GO_clusters(FBgn_list = split(rownames(mat), hc), 
               N_top = 10, 
               cex.balloons = 1.5, 
               padj_cutoff = 0.05)
# Network
par(mar= c(1,1,1,1), cex= 1)
vl_STRING_network(net,
                  cex.vertices = 2,
                  cex.vertices.labels = 0.5)
legend("topleft",
       bty= "n",
       fill= Cc,
       legend = paste0("cluster ", seq(max(hc))))
dev.off()

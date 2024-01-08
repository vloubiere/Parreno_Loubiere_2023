setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)

# Import metadata ----
dat <- data.table(FC_file= list.files("db/FC_tables/ATAC/", "^ATAC.*_vs_PH18.txt$", full.names = T))
dat[, cdition:= gsub("^ATAC_|_vs_PH18|.txt$", "", basename(FC_file))]
dat <- dat[, fread(FC_file), cdition]
dat[, c("seqnames", "start", "end"):= vl_toDTranges(ID)]

# Select peaks present in one of the 3 conditions and significantly affected
dat <- dat[ID %in% dat[abs(log2FoldChange)>1 & padj<0.001 & log10(baseMean)>1.25, ID]]

# Clip outliers ----
dat[, log2FoldChange_corr:= {
  lim <- quantile(log2FoldChange, c(0.05, 0.95), na.rm= T)
  log2FoldChange[log2FoldChange<lim[1]] <- lim[1]
  log2FoldChange[log2FoldChange>lim[2]] <- lim[2]
  log2FoldChange
}, cdition]

# Clustering using FC and padj----
layers <- dcast(dat, ID~cdition, value.var = c("log2FoldChange_corr", "padj"))
layers[, padj_PH29:= -log10(padj_PH29)]
layers[, padj_PHD11:= -log10(padj_PHD11)]
layers <- list(as.matrix(layers[,1:2], 1),
               as.matrix(layers[,c(1,3)], 1),
               as.matrix(layers[,c(1,4)], 1),
               as.matrix(layers[,c(1,5)], 1))
grid <- somgrid(1, 
                3, 
                "hexagonal", 
                toroidal= T)
set.seed(15)
som <- supersom(data = layers, 
                grid= grid)
som$unit.classif <- c(3,1,2)[som$unit.classif]

# SAVE ----
saveRDS(som,
        "Rdata/clustering_ATAC.rds")

# Plot ----
mat <- dcast(dat, ID~cdition, value.var = "log2FoldChange")
mat <- as.matrix(mat, 1)
hm <- vl_heatmap(mat, 
                 row.clusters = som$unit.classif,
                 cluster.rows= F,
                 cluster.cols= F, 
                 breaks = c(-3, 0, 3), 
                 col= c("cornflowerblue", "white", "tomato"),
                 show.rownames = F, 
                 legend.title = "Fold Change (log2)", 
                 show.col.clusters = F)

par(mar= c(5.1, 4.1, 4.1, 9.1),
    las= 2)
plot(hm)
.t <- rev(table(som$unit.classif))
axis(2,
     cumsum(.t)-.t/2,
     .t,
     lty= 0)

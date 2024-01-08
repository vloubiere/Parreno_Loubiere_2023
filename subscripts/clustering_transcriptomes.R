setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)

# Import metadata ----
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_noGFP" & FC_file!="NA"]
dat <- meta[, fread(FC_file), .(DESeq2_object, cdition, FC_file)]
# any signif foldChange
diffGenes <- dat[diff %in% c("up", "down"), FBgn]
dat <- dat[FBgn %in% diffGenes] 
# Not affected in ph18 (no FC cutoff! and remove ph18)
diffGenesPH18 <- dat[padj<0.05 & cdition=="PH18", FBgn]
dat <- dat[!(FBgn %in% diffGenesPH18)]
dat <- dat[cdition!="PH18"]
dat <- dat[, cdition:= factor(cdition, c("PHD11", "PHD9", "PH29"))]

# Clip outliers ----
dat[, corr_log2FoldChange:= {
  lim <- quantile(log2FoldChange, c(0.05, 0.95), na.rm= T)
  log2FoldChange[log2FoldChange<lim[1]] <- lim[1]
  log2FoldChange[log2FoldChange>lim[2]] <- lim[2]
  scale(log2FoldChange)
}, FC_file]
dat[padj>0.05, c("log2FoldChange", 
                 "corr_log2FoldChange"):= .(NA, 0)]

# Clustering ----
layers <- as.matrix(dcast(dat, FBgn~cdition, value.var = "corr_log2FoldChange"), 1)
layers <- list(constant= layers[, "PH29", drop= F],
               transient= layers[, c("PHD11", "PHD9")])
grid <- somgrid(3, 
                2, 
                "hexagonal", 
                toroidal= T)
init <- lapply(layers, function(x)
{
  set.seed(10)
  x <- x[sample(nrow(x), grid$xdim*grid$ydim), , drop= F]
  return(x)
})
som <- supersom(data = layers, 
                grid= grid,
                init = init,
                user.weights= c(1,1), 
                maxNA.fraction = 1)
som$unit.classif <- c(5, 4, 2, 6, 1, 3)[som$unit.classif]

mat <- as.matrix(dcast(dat, FBgn~cdition, value.var = "log2FoldChange"), 1)
vl_heatmap(mat[, c("PH29", "PHD9", "PHD11")], 
           row.clusters = som$unit.classif,
           cluster.rows= F,
           cluster.cols= F, 
           breaks = c(-4, 0, 4), 
           show.rownames = F, 
           legend.title = "Fold Change (log2)", 
           show.col.clusters = F)

# SAVE
saveRDS(som,
        "Rdata/clustering_RNA.rds")


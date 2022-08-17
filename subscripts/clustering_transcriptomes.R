setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)

##########################################################
# Import metadata
##########################################################
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA" & FC_file!="NA"]
dat <- meta[, fread(FC_file), .(DESeq2_object, cdition, FC_file)]
# any signif foldChange
diffGenes <- dat[diff %in% c("up", "down"), FBgn]
dat <- dat[FBgn %in% diffGenes] 
# Not affected in ph18 (no FC cutoff! and remove ph18)
diffGenesPH18 <- dat[padj<0.05 & cdition=="RNA_PH18", FBgn]
dat <- dat[!(FBgn %in% diffGenesPH18)]
dat <- dat[cdition!="RNA_PH18"]
dat <- dat[, cdition:= factor(cdition, c("RNA_PHD11", "RNA_PHD9", "RNA_PH29"))]

# Clip outliers
dat[, corr_log2FoldChange:= {
  lim <- quantile(log2FoldChange, c(0.05, 0.95), na.rm= T)
  log2FoldChange[log2FoldChange<lim[1]] <- lim[1]
  log2FoldChange[log2FoldChange>lim[2]] <- lim[2]
  scale(log2FoldChange)
}, FC_file]
dat[padj>0.05, c("log2FoldChange", 
                 "corr_log2FoldChange"):= .(NA, 0)]

##########################################################
# Clustering
##########################################################
layers <- as.matrix(dcast(dat, FBgn~cdition, value.var = "corr_log2FoldChange"), 1)
layers <- list(constant= layers[, "RNA_PH29", drop= F],
               transient= layers[, c("RNA_PHD11", "RNA_PHD9")])
grid <- somgrid(2, 
                3, 
                "hexagonal", 
                toroidal= T)
init <- lapply(layers, function(x)
{
  set.seed(8)
  x <- x[sample(nrow(x), grid$xdim*grid$ydim), , drop= F]
  return(x)
})
som <- supersom(data = layers, 
                grid= grid,
                init = init,
                user.weights= c(1,1), 
                maxNA.fraction = 1)

# SAVE
saveRDS(som,
        "Rdata/clustering_RNA.rds")
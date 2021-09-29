require(vlfunctions)
require(kohonen)

# Import conditions of interest
dat <- readRDS("Rdata/final_FC_table.rds")

# Select genes differentially expressed in at least one of the 6 conditions
dat[is.na(padj), padj:= 1]
dat <- dat[, .(!any(padj[cdition %in% c("EZ18_TS", "PH18_TS")]<0.05) & 
                 any(padj[!cdition %in% c("EZ18_TS", "PH18_TS")]<0.01),
               log2FoldChange, 
               cdition) , FBgn]
dat <- dat[(V1) & !cdition %in% c("EZ18_TS", "PH18_TS")]

# Clip extremes
vl_boxplot(dat, log2FoldChange~cdition, outline = T)
abline(h= c(-8, 8))
dat[, clipped_log2FoldChange:= log2FoldChange]
dat[log2FoldChange> 8, clipped_log2FoldChange:= 8]
dat[log2FoldChange< -8, clipped_log2FoldChange:= -8]

# Clustering
grid <- somgrid(3, 3, "hexagonal", toroidal= T)
train <- as.matrix(dcast(dat, 
                         FBgn~cdition, 
                         value.var = "clipped_log2FoldChange"), 1)
set.seed(5)
som <- supersom(train, 
                grid= grid, 
                maxNA.fraction = 0.1, 
                user.weights = c(1, 10))

# Add cluster columns to the data
som_cl <- as.data.table(som$data, keep.rownames = "FBgn")
som_cl[, cl:= som$unit.classif]
dat[som_cl, cl:= i.cl, on= "FBgn"]
dat <- na.omit(dat)
setorderv(dat, c("cl", "FBgn"))

# Add genes symbols
symbols <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
dat[symbols, symbol:= i.gene_symbol, on= "FBgn==gene_id"]

# Save tables and som
saveRDS(som, "Rdata/som_clustering_transcriptomes.rds")
saveRDS(dat, "Rdata/final_clustering_transcriptomes.rds")

mat <- dcast(dat, 
             cl+FBgn~cdition, 
             value.var = "clipped_log2FoldChange")
vl_heatmap(as.matrix(mat[, -1], 1), 
           cluster_rows = F, 
           show_rownames = F, 
           breaks = c(-5, 0, 5))

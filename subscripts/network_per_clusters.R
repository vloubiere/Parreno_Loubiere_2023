setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(cl)]

# Compute gene network
cols <- c("log2FoldChange_PH29", "log2FoldChange_PHD9", "log2FoldChange_PHD11")
sizes <- apply(dat[, ..cols], 1, function(x) max(abs(x), na.rm=T))
sizes[sizes>quantile(sizes, 0.90)] <- quantile(sizes, 0.90)
sizes <- log2(sizes+2)
dat[, sizes:= sizes]

# Clusters
g1 <- split(dat, dat$cl)
names(g1) <- paste0("cluster ", names(g1))
net1 <- lapply(g1, function(x) {
  x[, vl_STRING_interaction(symbols = symbol,
                            species = "Dm",
                            col= ifelse(PRC1_bound, "tomato", "cornflowerblue"),
                            size = sizes*2.5,
                            cex.label = sizes/4,
                            plot= F)]
})

#----------------------------------------------------#
# PLOT
#----------------------------------------------------#
pdf("pdf/cluster_network_per_cluster.pdf", height = 5, width = 5)
par(mar= c(2,2,2,2))
lapply(seq(net1), function(i)
{
  plot(net1[[i]],
       score_cutoff= 600, 
       top_N= 200)
  legend("topleft",
         bty= "n", 
         fill= c("tomato", "cornflowerblue"),
         legend= c("PRC1 +",
                   "PRC1-"))
  title(main= names(g1)[i])
})
dev.off()

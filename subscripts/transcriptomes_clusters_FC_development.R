setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

##########################
# Import data
##########################
# Dev transcriptomes
dev <- list.files("db/FC_tables/RNA/", 
                  pattern = "vs_RNA_72hED.txt$", 
                  full.names = T)
names(dev) <- gsub("RNA_development_Public_|RNA_|.txt$", "", basename(dev))
dev <- rbindlist(lapply(dev, fread), 
                 idcol = "variable")
# Add clusters
dat <- readRDS("Rdata/clustering_RNA.rds")$data
dev[dat, c("cl", "PRC1"):= .(i.cl, !is.na(i.PRC1_cluster)), on= "FBgn"]
dev <- dev[!is.na(cl)]

##########################
# Plot
##########################
pdf("pdf/Figures/Cluster_PRC1_bound_unbound_development_RNA.pdf",
    width= 9,
    height= 6)
par(las= 1,
    mar= c(2,3.5,4,0.5),
    mgp= c(2.25, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4))
dev[, {
  vl_boxplot(log2FoldChange~PRC1+cl,
             .SD,
             boxcol= c("lightgrey", "tomato"),
             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
             xaxt= "n",
             ylab= "log2FoldChange",
             ylim= c(-10, 10))
  title(main= variable, 
        line=3)
  abline(h= 0, lty= "11")
  axis(1, 
       at = seq(1.5, 11.5, 2), 
       labels = paste0("cl ", 1:6))
  legend(6,
         par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
         legend= c("noPRC1", "PRC1 bound"),
         fill= c("lightgrey", "tomato"),
         bty= "n",
         xpd=T)
  print("done")
}, variable]
dev.off()
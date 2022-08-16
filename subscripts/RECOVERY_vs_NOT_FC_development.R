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
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
dev[dat, c("PHD9_RECOVERY", "PHD11_RECOVERY"):= .(i.PHD9_RECOVERY, i.PHD11_RECOVERY), on= "FBgn"]
dev <- melt(dev, 
            id.vars = c("FBgn", "variable", "log2FoldChange"),
            measure.vars = c("PHD9_RECOVERY", "PHD11_RECOVERY"), 
            variable.name = "recovery")
dev <- na.omit(dev)
dev[, fisher.test(value, log2FoldChange>1, alternative = "greater")[c("estimate", "p.value")], .(variable, recovery)]
dev[, fisher.test(value, log2FoldChange<(-1), alternative = "greater")[c("estimate", "p.value")], .(variable, recovery)]

##########################
# Plot
##########################
pdf("pdf/Figures/dev_FC_revert_vs_not.pdf",
    width= 5,
    height= 6)
par(las= 2,
    mar= c(2,3.5,4,0.5),
    mgp= c(2.25, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4),
    cex.main= 0.5)
dev[, {
  vl_boxplot(log2FoldChange~value,
             .SD,
             boxcol= c("lightgrey", "tomato"),
             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
             ylab= "log2FoldChange",
             ylim= c(-5, 5))
  title(main= paste(recovery, variable), 
        line=3)
  abline(h= 0, lty= "11")
  legend(6,
         par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
         legend= c("noPRC1", "PRC1 bound"),
         fill= c("lightgrey", "tomato"),
         bty= "n",
         xpd=T)
  print("done")
}, .(variable, recovery)]
dev.off()
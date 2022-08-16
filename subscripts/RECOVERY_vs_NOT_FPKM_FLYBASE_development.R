setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

##########################
# Import data
##########################
# FLYBASE tissue transcriptomes
dev <- fread("/mnt/d/_R_data/genomes/dm6/gene_rpkm_report_fb_2021_02.tsv", skip= 5)
dev <- dev[Parent_library_name=="modENCODE_mRNA-Seq_development"]
dev <- dev[RNASource_name %in% c("mE_mRNA_em14-16hr",
                                 "mE_mRNA_em16-18hr",
                                 "mE_mRNA_em18-20hr",
                                 "mE_mRNA_em20-22hr",
                                 "mE_mRNA_em22-24hr",
                                 "mE_mRNA_L1",
                                 "mE_mRNA_L2",
                                 "mE_mRNA_L3_12hr",
                                 "mE_mRNA_WPP")]
dev[, RNASource_name:= factor(RNASource_name, c("mE_mRNA_em14-16hr",
                                                "mE_mRNA_em16-18hr",
                                                "mE_mRNA_em18-20hr",
                                                "mE_mRNA_em20-22hr",
                                                "mE_mRNA_em22-24hr",
                                                "mE_mRNA_L1",
                                                "mE_mRNA_L2",
                                                "mE_mRNA_L3_12hr",
                                                "mE_mRNA_WPP"))]
# Add clusters
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
dev[dat, c("PHD9_RECOVERY", "PHD11_RECOVERY"):= .(i.PHD9_RECOVERY, i.PHD11_RECOVERY), on= "FBgn#==FBgn"]
dev <- melt(dev, 
            id.vars = c("FBgn#", "RNASource_name", "RPKM_value"),
            measure.vars = c("PHD9_RECOVERY", "PHD11_RECOVERY"), 
            variable.name = "recovery")
dev <- na.omit(dev)

##########################
# Plot
##########################
pdf("pdf/Figures/FLYBASE_dev_FPKM_revert_vs_not.pdf",
    width= 10,
    height= 6)
par(las= 2,
    mar= c(2,3.5,4,0.5),
    mgp= c(2.25, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4),
    cex.main= 0.5)
dev[, {
  vl_boxplot(RPKM_value~value+RNASource_name,
             .SD,
             boxcol= c("lightgrey", "tomato"),
             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
             ylab= "log2FoldChange")
  title(main= paste(recovery), 
        line=3)
  abline(h= 0, lty= "11")
  legend(6,
         par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
         legend= c("noPRC1", "PRC1 bound"),
         fill= c("lightgrey", "tomato"),
         bty= "n",
         xpd=T)
  print("done")
}, recovery]
dev.off()
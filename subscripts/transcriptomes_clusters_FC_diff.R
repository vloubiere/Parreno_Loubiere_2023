setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import data
dat <- readRDS("Rdata/clustering_RNA.rds")$data
cols <- c("log2FoldChange_PHD11",
          "log2FoldChange_PHD9", 
          "log2FoldChange_PH29", 
          "PH18_FPKM",
          "H3K27me3_PH18",
          "H3K27me3_PHD11",
          "H3K27me3_PHD9",
          "H3K27me3_PH29",
          "H2AK118Ub_PH18",
          "H2AK118Ub_PHD11",
          "H2AK118Ub_PHD9",
          "H2AK118Ub_PH29",
          "H3K4me1_PH18",
          "H3K4me1_PHD11",
          "H3K4me1_PHD9",
          "H3K4me1_PH29",
          "H3K27Ac_PH18",
          "H3K27Ac_PHD11",
          "H3K27Ac_PHD9",
          "H3K27Ac_PH29",
          "H3K36me3_PH18",
          "H3K36me3_PHD11",
          "H3K36me3_PHD9",
          "H3K36me3_PH29",
          "PH_PH18",
          "PH_PHD11",
          "PH_PHD9",
          "PH_PH29",
          "PH_ED",
          "PC_ED",
          "SUZ12_ED",
          "H3K27me2_ED")
dat <- melt(dat, 
            id.vars = c("cl", "PRC1_cluster"), 
            measure.vars = cols)
dat[grepl("FPKM", variable), c("cdition", "ChIP"):= tstrsplit(variable, "_")]
dat[!grepl("FPKM", variable), c("ChIP", "cdition"):= tstrsplit(variable, "_")]
dat[grepl("_ED$", variable), ChIP:= paste0(ChIP, "_SA")]
dat[, ylab:= switch(ChIP, 
                    "FPKM"= "FPKM",
                    "log2FoldChange"= "FoldChange (log2)",
                    "Enrichment"), ChIP]
dat[, PRC1:= !is.na(PRC1_cluster)]
dat[, c("ymin", "ymax"):= as.list(vl_boxplot(value~PRC1+cl,
                                             .SD,
                                             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
                                             plot= F)$ylim), variable]
dat[, ymin:= min(ymin), ChIP]
dat[, ymax:= max(ymax), ChIP]
dat[, variable:= factor(variable, cols)]
setkeyv(dat, "variable")

pdf("pdf/Figures/Cluster_PRC1_bound_unbound_RNA.pdf",
    width= 9,
    height= 6)
par(las= 1,
    mar= c(2,3.5,4,0.5),
    mgp= c(2.25, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4))
dat[, {
  vl_boxplot(value~PRC1+cl,
             .SD,
             boxcol= c("lightgrey", "tomato"),
             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
             xaxt= "n",
             ylim= c(ifelse(ymin<0, ymin, 0), ymax),
             ylab= ylab)
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
}, .(variable, ymin, ymax, ylab)]
dev.off()
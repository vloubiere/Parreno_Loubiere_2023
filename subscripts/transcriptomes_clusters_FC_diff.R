setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[is.na(cl), cl:= 0]
# Select columns to plot
cols <- c("log2FoldChange_PH18",
          "log2FoldChange_PH29", 
          "log2FoldChange_PHD9", 
          "log2FoldChange_PHD11",
          "FPKM_PH18",
          "FPKM_PH29",
          "FPKM_PHD9",
          "FPKM_PHD11",
          "H3K27me3_PH18_BODY",
          "H3K27me3_PH29_BODY",
          "H3K27me3_PHD9_BODY",
          "H3K27me3_PHD11_BODY",
          "H2AK118Ub_PH18_BODY",
          "H2AK118Ub_PH29_BODY",
          "H2AK118Ub_PHD9_BODY",
          "H2AK118Ub_PHD11_BODY",
          "H3K4me1_PH18_BODY",
          "H3K4me1_PH29_BODY",
          "H3K4me1_PHD9_BODY",
          "H3K4me1_PHD11_BODY",
          "H3K27Ac_PH18_BODY",
          "H3K27Ac_PH29_BODY",
          "H3K27Ac_PHD9_BODY",
          "H3K27Ac_PHD11_BODY",
          "H3K36me3_PH18_BODY",
          "H3K36me3_PH29_BODY",
          "H3K36me3_PHD9_BODY",
          "H3K36me3_PHD11_BODY",
          "PH_PH18_TSS",
          "PH_PH29_TSS",
          "PH_PHD9_TSS",
          "PH_PHD11_TSS",
          "PH_ED_TSS",
          "PC_ED_TSS",
          "SUZ12_ED_TSS")
# Melt data
.m <- melt(dat[, c("cl", "PRC1_bound", cols), with= F], 
           id.vars = c("cl", "PRC1_bound"), 
           measure.vars = cols)
# Transform fpkm to log
.m[grepl("FPKM", variable), value:= log2(value+0.001)]
# Make unique conditions to predict ylims
.m[, cdition:= tstrsplit(variable, "_", keep= 1)]
.m[grepl("ED", variable), cdition:= paste0(cdition, "_ChIP")]
.m[, ylab:= fcase(grepl("FPKM", cdition), "FPKM",
                  grepl("log2FoldChange", cdition), "FoldChange (log2)",
                  default= "Enrichment")]
.m[, c("ymin", "ymax"):= as.list(vl_boxplot(value~PRC1_bound+cl,
                                            .SD,
                                            compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
                                            plot= F)$ylim), variable]
.m[, ymin:= min(ymin), cdition]
.m[, ymax:= max(ymax), cdition]
.m[, variable:= factor(variable, cols)]
setkeyv(.m, "variable")

pdf("pdf/Figures/Cluster_PRC1_bound_unbound_RNA.pdf",
    width= 9,
    height= 6)
par(las= 2,
    mar= c(2,3.5,4,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4))
.m[, {
  vl_boxplot(value~PRC1_bound+cl,
             .SD,
             boxcol= c("lightgrey", "tomato"),
             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
             xaxt= "n",
             ylim= c(ifelse(ymin<0, ymin, 0), ymax),
             ylab= ylab)
  title(main= variable,
        line= 2.75)
  abline(h= 0,
         lty= "11")
  axis(1, 
       at = seq(1.5, 13.5, 2), 
       labels = paste0("cl ", 0:6))
  legend(par("usr")[1],
         par("usr")[4]+strheight("M")*4.5, 
         legend= c("PRC1-", "PRC1+"),
         fill= c("lightgrey", "tomato"),
         bty= "n",
         xpd=T)
  print("done")
}, .(variable, ymin, ymax, ylab)]
dev.off()
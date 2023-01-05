require(vlfunctions)

dat <- fread("Rdata/final_gene_features_table.txt")
pl <- melt(dat, 
           id.vars = "recovery", 
           measure.vars = c("log2FoldChange_PH18", "log2FoldChange_PH29", "log2FoldChange_PHD11",
                            "FPKM_PH18", "FPKM_PH29", "FPKM_PHD11",
                            "H3K27me3_PH18_body", "H3K27me3_PH29_body", "H3K27me3_PHD11_body",
                            "H3K27Ac_PH18_prom", "H3K27Ac_PH29_prom", "H3K27Ac_PHD11_prom",
                            "PH_PH18_prom", "PH_PH29_prom", "PH_PHD11_prom"), 
           variable.factor = F)
pl <- rbind(pl[!is.na(recovery)],
            pl[, recovery:= "all"])
pl[, recovery:= factor(recovery, c("all", "noRecovery", "Recovery"))]
pl[, class:= tstrsplit(variable, "_", keep= 1)]
pl[, class:= switch(class, 
                    "log2FoldChange"= "RNA-Seq fold change (log2)",
                    "FPKM"= "RNA-Seq FPKM",
                    "H3K27me3"= "H3K27me3 enrichment",
                    "H3K27Ac"= "H3K27Ac enrichment",
                    "PH"= "PH enrichment"), class]
pl[class=="PH enrichment", value:= value+1]

pdf("pdf/Figure_3_boxplots_FPKM_CUTNRUN_ChIP.pdf", 
    width = 5, 
    height = 4.5)
par(mar= c(10,5,3,10),
    las= 1,
    tcl= -0.2,
    mgp= c(1.5, 0.5, 0))
pl[, {
  stats <- if(class=="RNA-Seq FPKM")
    list(c(2,3), c(5,6), c(8,9), c(2,8), c(3,9)) else
      list(c(1,2), c(2,3), c(1,3), c(4,5), c(5,6), c(4,6), c(7,8), c(8,9), c(7,9))
  log <- ifelse(class=="PH enrichment", "y", "")
  vl_boxplot(value~recovery+variable, 
             .SD, 
             tilt.names= T, 
             compute_pval= stats, 
             col= c("lightgrey", "rosybrown1", "palegreen3"), 
             boxwex= 0.6,
             xaxt= "n",
             at= c(1,2,3,6,7,8,11,12,13),
             ylab= class,
             log= log)
  axis(1, at= c(2, 7, 12), labels = NA)
  vl_tilt_xaxis(x= c(2, 7, 12),
                y = grconvertY(grconvertY(0, "npc", "inch")-strheight("M", "inch"), "inch", "user"),
                labels = c("No ph-KD", "Constant ph-KD", "Transient ph-KD"),
                srt= 30)
  legend(par("usr")[2],
         par("usr")[4],
         fill= c("rosybrown1", "palegreen3", "lightgrey"),
         legend= c("Irreversible", "Reversible", "Unaffected genes"),
         bty= "n",
         xpd= T)
  print("")
}, class]
dev.off()
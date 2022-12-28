# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[!is.na(FC_peaks), fread(FC_peaks), .(ChIP, cdition, FC_peaks)]
dat[, col:= fcase(padj<0.05 & log2FoldChange>1, "tomato",
                  padj<0.05 & log2FoldChange<(-1), "cornflowerblue",
                  default= "lightgrey")]

pdf("pdf/cutnrun_MA_plots.pdf",
    height = 5.5)
par(mfrow=c(2,3),
    mgp= c(2,0.5,0),
    tcl= -0.2)
dat[, {
  # MA plot
  plot(log2(baseMean),
       log2FoldChange, 
       pch= 19,
       col= adjustcolor(col, 0.5),
       cex= 0.7,
       ylim= c(-3.5, 3.5),
       las= 1,
       main= paste(ChIP, cdition))
  abline(h=0, lty= "11")
  legend("topright",
         c(paste0("Up (", sum(col=="tomato"), ")"),
           paste0("Unaffected (", sum(col=="lightgrey"), ")"),
           paste0("Down (", sum(col=="cornflowerblue"), ")")),
         col= c("tomato", "lightgrey", "cornflowerblue"),
         pch= 19,
         bty= "n",
         cex= 0.5)
  print("")
}, .(ChIP, cdition)]
dev.off()

png("pdf/cutnrun_MA_plots.png",
    res = 600,
    width= 4000,
    height = 3140*3)
par(mfrow=c(6,3),
    mgp= c(2,0.5,0),
    tcl= -0.2)
dat[, {
  # MA plot
  plot(log2(baseMean),
       log2FoldChange, 
       pch= 19,
       col= adjustcolor(col, 0.5),
       cex= 0.7,
       ylim= c(-3.5, 3.5),
       las= 1,
       main= paste(ChIP, cdition))
  legend("topright",
         c(paste0("Up (", sum(col=="tomato"), ")"),
           paste0("Unaffected (", sum(col=="lightgrey"), ")"),
           paste0("Down (", sum(col=="cornflowerblue"), ")")),
         col= c("tomato", "lightgrey", "cornflowerblue"),
         pch= 19,
         bty= "n",
         cex= 0.5)
  abline(h=0, lty= "11")
  print("")
}, .(ChIP, cdition)]
dev.off()

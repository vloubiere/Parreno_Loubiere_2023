setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(ChIP, cdition, FC_file)]

# Check if overlaps with reversed mark
dat[ChIP %in% c("H3K27Ac", "H3K27me3"), 
    c("peaks_file", "reversed"):= 
      {
        file <- list.files("db/FC_tables/cutnrun/", 
                           switch(ChIP, 
                                  "H3K27me3"= paste0("H3K27Ac_", cdition, "_vs_PH18.txt"),
                                  "H3K27Ac"= paste0("H3K27me3_", cdition, "_vs_PH18.txt")),
                           full.names = T)
        if(ChIP=="H3K27me3")
          peaks <- fread(file)[padj<0.05 & log2FoldChange>log2(1.5)] 
        if(ChIP=="H3K27Ac")
          peaks <- fread(file)[padj<0.05 & log2FoldChange<(-log2(1.5))]
        .(file, 
          peaks[.SD, .N>0, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1)
      }, .(ChIP, cdition)]


pdf("pdf/cutnrun/MA_plots.pdf", 
    height = 5.5)
par(mfrow=c(2,3))
dat[, {
  # MA plot
  Cc <- fcase(padj<0.05 & log2FoldChange>log2(1.5), "tomato",
              padj<0.05 & log2FoldChange<(-log2(1.5)), "cornflowerblue",
              default= "lightgrey")
  plot(log2(baseMean),
       log2FoldChange, 
       pch= 19,
       col= adjustcolor(Cc, 0.5),
       cex= 0.7,
       ylim= c(-3.5, 3.5),
       las= 1,
       main= paste(ChIP, cdition))
  abline(h=0, lty= "11")
  # Add reciprocal changes K27
  if(ChIP %in% c("H3K27me3", "H3K27Ac"))
  {
    points(log2(baseMean)[reversed],
           log2FoldChange[reversed],
           cex= 0.5)
    legend(switch(ChIP, 
                  "H3K27me3"= "topleft",
                  "H3K27Ac"= "bottomleft"),
           legend = c(paste0("Up ", sum(Cc=="tomato")),
                      paste0("Stable ", sum(Cc=="lightgrey")),
                      paste0("Down ", sum(Cc=="cornflowerblue")),
                      paste0(switch(ChIP,
                                    "H3K27Ac" = "H3K27me3 Loss ",
                                    "H3K27me3" = "H3K27Ac Gain "), sum(reversed))),
           pch= c(19,19,19,1),
           col= c(adjustcolor(c("tomato", "lightgrey", "cornflowerblue"), 0.5), "black"),
           text.col= c("tomato", "grey20", "cornflowerblue", "black"),
           bty= "n",
           cex= 0.5)
  }
  print("")
}, .(ChIP, cdition)]
dev.off()

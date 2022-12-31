# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(png)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[grepl("H2AK118Ub", ChIP) 
            & !is.na(FC_peaks)
            & cdition %in% c("PH29", "PHD11"), fread(FC_peaks), .(ChIP, cdition, FC_peaks)]
dat[, diff:= fcase(padj<0.05 & log2FoldChange>1, "up",
                   padj<0.05 & log2FoldChange<(-1), "down",
                   default= "unaffected")]
dat[, col:= switch(diff, "up"= "tomato", "down"= "cornflowerblue", "unaffected"= "grey"), diff]
dat[, cdition:= switch(cdition, "PH29"= "Constant ph-KD", "PHD11"= "Transient ph-KD"), cdition]

pdf("pdf/Extended_data_5_MA_plots_K118Ub.pdf",
    height = 6,
    width = 5.5)
par(mfrow=c(2,2),
    mgp= c(2,0.5,0),
    tcl= -0.2,
    las= 1)
dat[, {
  tmp <- tempfile(fileext = "png")
  png(tmp,
      res = 600,
      width= 4000,
      height = 4000,
      type="cairo")
  par(mar= c(0,0,0,0))
  plot.new()
  plot.window(xlim= range(log10(baseMean)), 
              ylim= c(-4, 4))
  points(log10(baseMean), 
         log2FoldChange,
         col= adjustcolor(col, 0.4),
         pch= 16,
         cex= 2)
  xlim <- par("usr")[1:2]
  ylim <- par("usr")[3:4]
  dev.off()
  
  pic <- readPNG(tmp)
  plot(NA,
       ylab= "log2FoldChange", 
       xlab= "baseMean",
       xlim= 10^xlim, 
       ylim= ylim,
       log= "x",
       frame= F)
  rasterImage(pic, 
              10^xlim[1], 
              ylim[1],
              10^xlim[2], 
              ylim[2])
  title(main= paste(ChIP, cdition), 
        font.main = 1)
  leg <- c(paste(sum(diff=="up"), "Up"), 
           paste(sum(diff=="unaffected"), "Unaffected"), 
           paste(sum(diff=="down"), "Down"))
  legend("topright", 
         leg, 
         bty= "n", 
         text.col = c("red", "grey", "blue"), 
         cex = 0.5)
  print("DONE")
}, .(ChIP, cdition)]
dev.off()

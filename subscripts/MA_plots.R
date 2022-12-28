# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(png)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file) & system=="noGFP", fread(FC_file), .(FC_file, cdition)]
dat[, cdition:= switch(cdition, 
                       "PH18"= "no ph-KD",
                       "PH29"= "Constant ph-KD",
                       "PHD9"= "Transient ph-KD d9",
                       "PHD11"= "Transient ph-KD d11"), cdition]
dat[, cdition:= factor(cdition, 
                        levels= c("no ph-KD", "Constant ph-KD", "Transient ph-KD d9", "Transient ph-KD d11"))]

setorderv(dat, "cdition")
dat[, col:= switch(diff, "up"= "tomato", "down"= "cornflowerblue", "unaffected"= "lightgrey"), diff]
dat <- dat[order(col=="lightgrey", decreasing= T)]

pdf("pdf/Figure_2_MA_plots_transcriptomes.pdf", 10, 2.75)
par(mfrow=c(1,4),
    tcl= -0.2,
    mgp= c(2,0.5,0),
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
              ylim= c(-10, 10))
  points(log10(baseMean), 
         ifelse(abs(log2FoldChange)>10, sign(log2FoldChange)*10, log2FoldChange),
         col= adjustcolor(col, 0.4),
         pch= ifelse(abs(log2FoldChange)>10, 17, 16),
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
  title(main= cdition, font.main = 1)
  leg <- c(paste(sum(diff=="up"), "Up"), 
           paste(sum(diff=="down"), "Down"), 
           "(p<0.05, |log2FC|>1)")
  legend("topright", 
         leg, 
         bty= "n", 
         text.col = c("red", "blue", "black"), 
         cex = 0.5)
  print("DONE")
}, .(cdition, FC_file)]
dev.off()

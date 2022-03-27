setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
options(scipen = 999)
require(kohonen)
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm6)

meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[!is.na(FC_file_ChIP_peaks), fread(.BY[[1]], colClasses = c("character", rep("numeric", 5))), .(FC_file_ChIP_peaks, ChIP, cdition)]

pdf("pdf/cutnrun/volcano_plots_K27_changes.pdf", 
    height = 5.25)
par(mfrow= c(2,3),
    las=1)
dat[, {
  xlim <- switch(ChIP, 
                 "H3K27me3"= c(-1.5,1.5),
                 "H3K27Ac"= c(-3,3))
  plot(log2FoldChange, 
       -log10(padj),
       main= paste(ChIP, cdition),
       xlim= xlim,
       ylim= c(0,40),
       pch= 16,
       col= adjustcolor(ifelse(padj<0.001, "tomato", "grey"), 0.4),
       cex= 0.5)
  print("DONE")
}, keyby= .(ChIP, cdition)]
dev.off()
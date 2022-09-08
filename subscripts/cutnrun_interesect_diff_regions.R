setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[!is.na(FC_peaks), fread(FC_peaks), .(ChIP, cdition, FC_peaks)]
dat[padj<0.05 & abs(log2FoldChange)>1, change:= paste0(ChIP, "_", cdition, "_",  ifelse(log2FoldChange>0, "up", "down")), .(cdition, ChIP)]
dat <- dat[grepl("H3K27me3.*down|H3K27Ac.*up|H3K4me1.*up", change)]
dat <- vl_importBed(dat[, c("seqnames", "start", "end"):= tstrsplit(ID, "_|:|-", keep= 2:4)])

.m <- vl_collapseBed(dat)
dat[.m, idx:= i.idx, on= c("seqnames", "start<=end", "end>=start")]
categ <- split(dat$idx, dat$change)

vl_upset_plot(pl[grep("PH29", names(categ), value = T)])
vl_upset_plot(pl[grep("H3K27me3", names(categ), value = T)])
vl_upset_plot(pl[grep("H3K27Ac", names(categ), value = T)])
vl_upset_plot(pl[grep("H3K4me1", names(categ), value = T)])


pdf("pdf/cutnrun/MA_plots.pdf",
    height = 5.5)
par(mfrow=c(2,3))
dat[, {
  # MA plot
  Cc <- fcase(padj<0.05 & log2FoldChange>1, "tomato",
              padj<0.05 & log2FoldChange<(-1), "cornflowerblue",
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
  print("")
}, .(ChIP, cdition)]
dev.off()

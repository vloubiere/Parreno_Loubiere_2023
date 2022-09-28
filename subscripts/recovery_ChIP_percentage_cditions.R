setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")[!is.na(recovery)]
cols <- c("FBgn", "recovery", grep("_body$|_prom$|_TTS$", names(dat), value = T))
dat <- dat[, ..cols]

# Compute percentage
HTM <- melt(dat, id.vars= c("FBgn", "recovery"))
HTM[, c("ChIP", "variable", "rep"):= tstrsplit(variable, "_")]
HTM <- merge(HTM[variable=="PH18"],
             HTM[variable!="PH18"], by= c("FBgn", "rep", "ChIP", "recovery"))
HTM[, perc:= value.y/value.x*100]
HTM <- HTM[, .(perc= mean(perc)), .(FBgn, ChIP, recovery, cdition= variable.y)]
HTM[, cdition:= factor(cdition, c("PH29", "PHD9", "PHD11"))]

pdf("pdf/recovery_ChIP_percentage_cditions.pdf", 
    width= 8,
    height= 2.5)
par(las= 2,
    mar= c(4,3.5,4,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(1,6))
HTM[, {
  vl_boxplot(perc~recovery+cdition,
             col= adjustcolor(c("palegreen3", "rosybrown1"), 0.7),
             compute_pval= list(c(1,2), c(3,4), c(5,6)),
             ylab= "% PH18 signal",
             xaxt= "n",
             at= rep(seq(1,5,2), each=2)+c(0.2,0.8))
  legend(par("usr")[1]-strwidth("M"),
         par("usr")[4]+strheight("M")*5,
         fill= c("palegreen3", "rosybrown1"),
         legend= c("Recovery", "No recovery"),
         bty= "n",
         xpd= T)
  abline(h= 100, lty= 2)
  axis(1, 
       seq(1.5, 5.5, 2),
       levels(cdition))
  title(main= ChIP, line= 3.1)
}, ChIP]
dev.off()

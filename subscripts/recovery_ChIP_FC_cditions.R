setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# FC
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
meta <- melt(meta, id.vars = c("ChIP", "cdition"), measure.vars = patterns("^FC_"))
meta <- meta[grepl("_body|_prom|_TTS", variable)]
meta <- na.omit(unique(meta))
FC <- meta[, fread(value), .(ChIP, cdition)]

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")[!is.na(recovery), .(FBgn, recovery)]
dat <- FC[dat, on= "ID==FBgn"]
dat[, cdition:= factor(cdition, c("PH29", "PHD9", "PHD11"))]
dat <- dat[!is.na(log2FoldChange)]

pdf("pdf/recovery_ChIP_FC_cditions.pdf", 
    width= 8,
    height= 2.5)
par(las= 2,
    mar= c(4,3.5,4,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(1,6))
dat[, {
  vl_boxplot(log2FoldChange~recovery+cdition,
             col= adjustcolor(c("palegreen3", "rosybrown1"), 0.7),
             compute_pval= list(c(1,2), c(3,4), c(5,6)),
             ylab= "FoldChange vs PH18 (log2)",
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
  abline(h= 0, lty= 2)
  title(main= ChIP, 
        line= 3.1)
}, ChIP]
dev.off()

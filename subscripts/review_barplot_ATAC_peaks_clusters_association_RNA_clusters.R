setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import peaks and TSSs ----
peaks <- readRDS("Rdata/final_ATAC_table.rds")
peaks <- vl_resizeBed(peaks, "center", 0, 0)
genes <- readRDS("Rdata/final_gene_features_table.rds")[!is.na(cl)]
genes[, cl:= factor(cl, c("Unaffected", "Transient-specific", "Irreversible", "Reversible", "Down 1", "Down 2", "Down 3"))]
TSS <- vl_resizeBed(genes, "start", 0, 0)

# Assign peaks to closest genes ----
dat <- vl_closestBed(peaks, TSS)
dat <- dat[abs(dist)<25e3]

# Fisher tests ----
cmb <- CJ(dat$cl,
          dat$cl.b,
          unique = T)
cmb <- cmb[V1!="Unaffected" & V2!="Unaffected"]
.f <- cmb[, {
  fisher.test(table(dat$cl==V1,
                    dat$cl.b==V2),
              alternative = "greater")[c("estimate", "p.value")]
}, .(V1, V2)]
.f[, padj:= p.adjust(p.value, "fdr")]
.f[, V2:= gsub("Cluster ", "", V2)]

# Plot ----
Cc <- colorRampPalette(c("grey10", "grey80"))(length(unique(.f$V2)))

pdf("pdf/review_ATAC_clusters_OR_RNA_cl_barplot.pdf", 4.1, 2)
vl_par(mfrow= c(3,1),
       mai= c(0,0,0.05,0),
       omi= c(.6, 1.6, .2, 1.8),
       cex.axis= 5/12,
       cex= 1,
       mgp= c(.5, .15, 0),
       tcl= -0.05)
.f[, {
  # Barplots
  bar <- barplot(estimate,
                 border= NA,
                 lwd= .5,
                 xaxt= "n",
                 yaxt= "n",
                 col= adjustcolor(Cc, .6),
                 width= 1)
  # pvalues
  vl_plot_pval_text(bar,
                    estimate,
                    padj,
                    stars.only = T,
                    xpd= NA,
                    cex= .5,
                    offset = 0,
                    show.NS = F)
  # Legends
  axis(4, lwd= .5)
  text(par("usr")[1],
       mean(par("usr")[c(3,4)]),
       V1,
       xpd= NA,
       pos= 2,
       cex= 9/12,
       offset= 0)
  segments(bar[1]-.5,
           1,
           bar[length(bar)]+.5,
           1,
           lty= "11",
           lwd= .5,
           lend= 3)
  if(.GRP==2)
    text(par("usr")[2]+strwidth("M")*1.25,
         mean(par("usr")[3:4]),
         "Odd ratio",
         xpd= NA,
         srt= -90)
  if(.GRP==.NGRP)
  {
    par(cex.axis= 9/12)
    vl_tilt_xaxis(bar, labels= V2)
  }
}, V1]
dev.off()
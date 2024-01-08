setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import classified K27 domains ----
K27 <- readRDS("Rdata/K27_domains_classif_reversible_irreversible.rds")

# Import log2FC ----
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
meta <- meta[ChIP %in% c("H3K27me3", "H3K27Ac") & cdition %in% c("PH29", "PHD11")]
FC <- meta[, fread(unique(FC_peaks)), .(ChIP, cdition)]
FC[, c("seqnames", "start", "end"):= tstrsplit(ID, "_|:|-", keep= 2:4, type.convert = T)]

# Compute overlaps ----
FC[, All:= vl_covBed(FC, K27)>0] # All log2FoldChanges used as control group
FC[, Irreversible:= vl_covBed(FC, K27[irreversible>0])>0]
FC[, Reversible:= vl_covBed(FC, K27[reversible>0])>0]

# melt ----
pl <- melt(FC,
           id.vars= c("ChIP", "ID", "cdition", "log2FoldChange"),
           measure.vars= c("All", "Irreversible", "Reversible"))
pl <- pl[(value)]
pl[, cdition:= switch(cdition,
                      "PH29"= "Constant\nph-KD",
                      "PHD11"= "Transient\nph-KD"), cdition]
pl[, cdition:= factor(cdition, c("Constant\nph-KD", "Transient\nph-KD"))]

# Plots ----
Cc <- c("lightgrey", "rosybrown1", "palegreen3")
Cc <- adjustcolor(Cc, .6)
at <- c(1:3, 5:7)

pdf("pdf/review_K27_domains_H3K27me3_H3K27Ac_FC.pdf",
    3,
    3)
vl_par(lwd= 0.5)
pl[, {
  # Boxplot
  box <- vl_boxplot(log2FoldChange~variable+cdition,
                    at= at,
                    col= Cc,
                    compute.pval= list(c(1,2), c(1,3), c(2,3), c(4,5), c(4,6), c(5,6)),
                    tilt.names= T,
                    ylab= "log2FoldChange",
                    main= ChIP,
                    xaxt= "n",
                    pval.cex= 5/12,
                    lwd= 0.5)
  abline(h= 0,
         lty= "11")
  axis(1,
       c(2,6),
       labels= levels(cdition),
       lwd= 0,
       cex= .8)
  # Legend (unique peaks per cartegory)
  unique(.SD[, .(variable, ID)])[, {
    leg <- paste0(levels(variable), " (", formatC(table(variable), big.mark = ","), ")")
    legend(par("usr")[2],
           par("usr")[4],
           legend= leg,
           fill= Cc,
           cex= 6/12,
           bty= "n",
           xpd= T)
  }]
  .SD
}, ChIP]
dev.off()
setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import FC ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[!is.na(cl)]
dat[, cl:= factor(cl, rev(levels(cl)))]
dat[, rescued:= diff_GFP_PH=="up" & diff_STAT92E_PH!="up"]
pl <- dat[, c(fisher.test(table(dat$cl==cl, dat$rescued),
                          alternative = "greater")[c("estimate", "p.value")],
              perc= sum(rescued, na.rm= T)/.N*100), keyby= cl]
pl[, padj:= p.adjust(p.value, "fdr")]

# Plot ----
Cc <- c("grey20", "grey50", "grey80")

pdf("pdf/review_cl_enrichment_Stat92E_rescued_genes.pdf",
    2.25,
    3)
vl_par(mgp= c(1.25,.35,0))
pl[, {
  bar <- barplot(estimate,
                 names.arg= cl,
                 beside= T,
                 horiz= T,
                 border= NA,
                 ylab= NA,
                 xaxt= "n",
                 xlab= "Rescued genes\nenrichment (odd ratio)")
  axis(1, padj = -1.75)
  abline(v= 1, lty= "11")
  perc <- paste0(round(perc, 1), "%")
  x <- estimate+strwidth("M")*1
  text(x,
       bar,
       perc,
       offset= 0.1,
       xpd= NA,
       cex= 6/12)
  vl_plot_pval_text(x,
                    bar,
                    padj,
                    stars.only = T,
                    show.NS = F,
                    xpd= NA,
                    pos= 3,
                    offset= .1)
}]
dev.off()


setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(glmnet)

# Import genes ----
dat <- readRDS("Rdata/final_ATAC_rescue_table.rds")
dat[, count:= ifelse(ZFH1_counts>quantile(ZFH1_counts, .99), quantile(ZFH1_counts, .99), ZFH1_counts)]

# Plot ----
pdf("pdf/review_ATAC_rescue_FC_vs_motif_counts.pdf", 2.5, 3)
vl_par(lwd= .75,
       mgp= c(.75, .25, 0))
dat[, {
  med <- vl_boxplot(log2FoldChange_zfh1_ph_vs_gfp_ph~count,
                    compute.pval= list(c(1, length(unique(count)))),
                    xaxt= "n",
                    xlab= "zfh1 motif counts",
                    ylab= "ATAC-Seq foldchange (log2)\nzfh1+ph-KD vs. gfp+ph-KD",
                    col= "lightgrey")
  axis(1,
       seq(unique(count)),
       sort(unique(count)),
       padj= -1.25)
  abline(h= med$stat[3,1], lty= "11")
  .SD
}]
dev.off()
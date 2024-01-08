setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(glmnet)

# Import genes ----
dat <- readRDS("Rdata/final_ATAC_table.rds")
dat <- dat[!is.na(log2FoldChange_PH29) & !is.na(log2FoldChange_PHD11)]

# motif counts
subCounts <- readRDS("db/counts/motifs/review_selected_ATAC_peaks_motif_counts_500.rds")
sel <- as.data.table(subCounts[, c("cad", "Stat92E", "zfh1")])
dat <- cbind(dat[, .(log2FoldChange_PH29, log2FoldChange_PHD11)], sel)
dat <- melt(dat,
            measure.vars = c("log2FoldChange_PH29", "log2FoldChange_PHD11"))
dat <- melt(dat,
            id.vars = c("variable", "value"))
setnames(dat, c("variable", "FC", "name", "count"))
dat[, cutoff:= quantile(count, .99), name]
dat[, count:= ifelse(count>cutoff, cutoff, count)]
dat[, Cc:= c("khaki3", "tomato", "lightblue")[.GRP], name]

# Plot
pdf("pdf/review_ATAC_FC_vs_motif_counts.pdf", width = 3.1, 4.6)
vl_par(mfrow= c(3,2),
       mai= c(.2,.05,.1,.1),
       omi= c(1,1,1,.5),
       lwd= .75,
       mgp= c(.5, .25, 0))
dat[, {
  med <- vl_boxplot(FC~count,
                    compute.pval= list(c(1, length(unique(count)))),
                    col= Cc,
                    xaxt= "n")
  axis(1,
       seq(unique(count)),
       sort(unique(count)),
       padj= -1.25)
  abline(h= med$stat[3,1], lty= "11")
  if(.GRP %% 2==1)
  {
    text(par("usr")[1],
         mean(par("usr")[3:4]),
         name,
         pos= 2,
         offset= 2,
         xpd= NA)
    title(ylab= "ATAC-Seq FC (log2)",
          line= .85,
          xpd= NA)
  }
  .SD
}, .(variable, name, Cc)]
dev.off()
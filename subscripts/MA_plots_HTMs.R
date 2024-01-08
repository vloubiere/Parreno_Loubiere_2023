setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(data.table)

# Import ----
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[grepl("^H3K27|^H2AK118Ub", ChIP) 
            & !is.na(FC_peaks)
            & cdition %in% c("PH29", "PHD11"), fread(FC_peaks), .(ChIP, cdition, FC_peaks)]
dat[, c("seqnames", "start", "end"):= tstrsplit(ID, "_|:|-", keep= 2:4, type.convert = T)]
dat[, diff:= fcase(padj<0.05 & log2FoldChange>1, "up",
                   padj<0.05 & log2FoldChange<(-1), "down",
                   default= "unaffected")]
dat[, col:= switch(diff, "up"= "tomato", "down"= "cornflowerblue", "unaffected"= "grey"), diff]
dat[, cdition:= switch(cdition, "PH29"= "Constant ph-KD", "PHD11"= "Transient ph-KD"), cdition]
dat[, diff:= factor(gsub("^(.)", "\\U\\1", diff, perl = TRUE),
                    c("Up", "Unaffected", "Down"))]
zfh1 <- readRDS("Rdata/final_gene_features_table.rds")[symbol=="zfh1"]
upd <- readRDS("Rdata/final_gene_features_table.rds")[symbol %in% c("upd1", "upd2", "upd3")]
dat[, zfh1:= vl_covBed(dat, zfh1)>0]
dat[, upd:= vl_covBed(dat, upd)>0]

# Plot ----
pdf("pdf/Figure_3_MA_plots_HTMs.pdf", 6, 6)
vl_par(mfrow=c(2,2))
dat[, {
  vl_rasterScatterplot(log10(baseMean),
                       log2FoldChange,
                       col= adjustcolor(col, .3),
                       xlab= "baseMean (log10)",
                       ylab= paste0(ChIP, " fold change (log2)"),
                       ylim =c(-4, 4),
                       pch= 16,
                       cex= .8,
                       main= cdition,
                       xaxt= "n")
  axis(1, padj= -1.5)
  .SD[(zfh1), points(log10(baseMean), log2FoldChange, cex= .5)]
  .SD[(upd), points(log10(baseMean), log2FoldChange, cex= .5, pch= 0)]
  legend("topright",
         paste0(levels(diff), " (", formatC(table(diff), big.mark = ","), ")"),
         bty= "n",
         text.col = c("red", "grey", "blue"),
         cex = 0.5,
         inset= c(0, -0.2),
         xpd= NA)
  legend("topleft",
         pch= c(1,0),
         legend= c("zfh1", "upd1-3"),
         bty= "n",
         cex = 0.5,
         inset= c(0, -0.2),
         xpd= NA)
  print("DONE")
}, .(ChIP, cdition)]
dev.off()
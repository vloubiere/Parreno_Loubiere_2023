setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data ----
SNP <- readRDS("Rdata/review_gDNA_final_table.rds")
# SNP <- SNP[allele_freq>0.2]
genes <- SNP[!is.na(symbol)]
genes <- genes[, .(symbol= unlist(tstrsplit(symbol, ","))), .(cdition, id)]
genes[, id:= paste0(id, "-", symbol)]
SV <- readRDS("Rdata/review_gDNA_SVs_final_table.rds")
CNV <- readRDS("Rdata/review_gDNA_CNVs_final_table.rds")
dat <- rbindlist(list(`SNV/InDels`= SNP,
                      Genes= genes,
                      SVs= SV,
                      CNVs= CNV),
                 idcol= "class",
                 fill= T)
dat <- unique(dat[, .(id, cdition, class)])

# plot
pdf("pdf/review_genetic_variants_overlaps.pdf", 5, 4*1.6)
vl_par(mfrow= c(4,1),
       mai= c(1, 1.2, .2, .2),
       mgp= c(0.4, 0.2, 0),
       cex.lab= .5,
       cex.axis= .3)
dat[, {
  .c <- vl_upset_plot(dat.list = split(id,
                                       cdition), 
                      grid.hex = .5,
                      cex.grid = .5,
                      cex.inter = .2,
                      intersection.cutoff = ceiling(.N/1000))
  title(main= class)
  sel <- .c$inter[, setdiff(names(.c$inter), c("N", "x")), with= F]
  perc <- sum(.c$inter[apply(sel, 1, sum)==1, N])
  perc <- perc/length(unique(id))
  text(grconvertX(0, "nfc", "user"),
       grconvertY(1, "nfc", "user")-strheight("M"),
       paste0(round(perc*100, 1), "% found in 1 sample"),
       pos= 4,
       xpd= NA,
       cex= 11/12)
  .SD
}, class]
dev.off()
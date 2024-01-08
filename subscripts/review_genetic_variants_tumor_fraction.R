setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data ----
# SNP
SNP <- readRDS("Rdata/review_gDNA_final_table.rds")
SNP <- SNP[allele_freq>0.2]
# Genes
genes <- SNP[!is.na(symbol) & !is.na(tumor_fraction)]
genes <- genes[, .(symbol= unlist(tstrsplit(symbol, ","))), .(cdition, tumor_fraction)]
setorderv(genes, "tumor_fraction", -1)
genes <- genes[, .(tumor_fraction= tumor_fraction[1]), .(id= symbol)]
# SVs
SV <- readRDS("Rdata/review_gDNA_SVs_final_table.rds")
SV <- SV[num_Reads>=5]
# CNVs
CNV <- readRDS("Rdata/review_gDNA_CNVs_final_table.rds")
CNV <- CNV[norm>1.5 | norm<2/3]
# Merge
dat <- rbindlist(list(`SNV/InDels`= SNP,
                      Genes= genes,
                      SVs= SV,
                      CNVs= CNV),
                 idcol= "class",
                 fill= T)
dat <- dat[!is.na(tumor_fraction), .(id, tumor_fraction, class)]
dat[, class:= factor(class, unique(class))]
dat <- unique(dat)

# Plot ----
pdf("pdf/review_SNP_InDel_tumor_fraction.pdf", 3, 3)
vl_par(mfrow= c(1,4),
       mai= c(1.1, 0, .9, .1),
       omi= c(0, .7, 0, .3),
       xpd= NA)
dat[, {
  # Barplot
  .t <- rev(table(tumor_fraction))
  bar <- barplot(.t,
          horiz= T,
          border= NA,
          yaxt= "n",
          xaxt= "n",
          col= c(rep("lightgrey", length(.t)-1), "tomato"))
  # Legend
  title(main= class, line= 1, font.main= 1, cex.main= 11/12)
  if(.GRP==1)
  {
    title(ylab= "Number of tumor samples", xpd= NA, line = 2)
    axis(2, bar, names(.t), lwd= 0)
  }
  # N single tumor
  Nsingle <- .t[length(.t)]
  perc <- round(Nsingle/sum(.t)*100, 1)
  axis(1, padj= -1.5)
  text(Nsingle/2,
       bar[length(bar)],
       paste0(perc, "%"),
       pos= 3,
       xpd= NA,
       cex= 7/12,
       offset= .75)
  text(Nsingle/2,
       bar[length(bar)],
       paste0("n= ", formatC(Nsingle, big.mark = ",")),
       pos= 3,
       xpd= NA,
       cex= 5/12,
       offset= .4)
  # N several tumors
  text(.t[-length(.t)],
       bar[-length(bar)],
       formatC(.t[-length(.t)], big.mark = ","),
       pos= 4,
       xpd= NA,
       cex= 4/12,
       offset= 0.2)
  .SD
}, class]
dev.off()
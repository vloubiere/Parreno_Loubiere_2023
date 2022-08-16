setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
dat <- dat[, .(FBgn, PHD9_RECOVERY, PHD11_RECOVERY)]
cl <- readRDS("Rdata/clustering_RNA.rds")$data
cl <- cl[, .(FBgn, cl, PH18_FPKM, PH29_FPKM, PHJ9_FPKM, PHJ11_FPKM)]
dat <- cl[dat, on= "FBgn"]

# melt data
pl <- melt(dat, measure.vars = patterns("RECOVERY$"))
pl <- pl[!is.na(value)]
pl[, cl:= factor(ifelse(is.na(cl), "NA", cl))]

#################################
# PLOT
#################################
pdf("pdf/Figures/compare_RECOVERY_vs_transcriptome_cl2&5.pdf",
    width= 10,
    height= 10)
par(mfrow= c(3,2))

# Plot counts across FC
pl[, {
  y <- switch(as.character(variable),
              "PHD9_RECOVERY"= as.matrix(SJ(PH18= PH18_FPKM, PH29= PH29_FPKM, PHJ9= PHJ9_FPKM)),
              "PHD11_RECOVERY"= as.matrix(SJ(PH18= PH18_FPKM, PH29= PH29_FPKM, PHJ11= PHJ11_FPKM)))
  pseudo <- round(min(y[y>0], na.rm = T), 3)
  y <- log2(y+pseudo)
  x <- matrix(seq(ncol(y)), 
              nrow = nrow(y), 
              ncol= ncol(y),  
              byrow = T)
  Cc <- ifelse(value, "limegreen", "red")
  plot(unlist(x), 
       unlist(y),
       col= Cc,
       pch= 16,
       xaxt= "n",
       xlab= NA,
       ylab= paste0("log2(FPKM+", pseudo, ")"))
  axis(1, 
       at= 1:3, 
       colnames(y), 
       las= 2)
  sapply(seq(nrow(y)), 
         function(i) 
           lines(x[i,], y[i,], col= Cc[i]))
  print("done")
}, variable]

# Transcriptome cluster fraction
pl[, {
  x <- table(cl, 
             useNA = "ifany")
  pie(x, labels = paste0("cl ", names(x), " (", x, ")"))
  title(main= paste0(variable, "= ", value))
}, .(variable, value)]

# Overlaps between D11 and D9
vl_upset_plot(split(pl$FBgn, pl[, .(variable, value)]))
dev.off()
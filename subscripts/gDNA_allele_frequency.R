setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!(PH18)]
pl <- unique(dat[, .(id, TUMOR_alt_freq_max, NORMAL_alt_freq_max, col, class, occurence)])

######################################################################
# Check link allele frequency and found in PH18
######################################################################
pdf("pdf/gDNA_allele_freq_scatterplots.pdf", 
    4, 
    4)
layout(matrix(1:2, ncol= 1), 
       heights = c(0.2, 1))
par(las= 1, 
    mgp= c(2, 0.5,0),
    tcl= -0.2)
pl[, {
  par(mar= c(0.2, 4.1, 1, 2.1))
  plot(density(NORMAL_alt_freq_max), 
       frame= F, 
       xaxt= "n", 
       main= paste0(occurence, " (SNP+InDel)"),
       yaxt= "n",
       ylab= NA)
  par(mar= c(3.1, 4.1, 0.2, 2.1))
  plot(NORMAL_alt_freq_max, TUMOR_alt_freq_max, 
       col= adjustcolor(col, 0.3),
       pch= 16,
       cex= 0.8)
  leg <- .SD[, .N, .(class, col)]
  legend("topright",
         legend= leg[, paste0(class, " (" , N, ")")],
         fill= leg$col,
         bty= "n")
  abline(0, 1)
}, occurence]
dev.off()

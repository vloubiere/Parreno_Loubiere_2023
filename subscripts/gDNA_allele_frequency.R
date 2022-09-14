setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- fread("Rdata/gDNA_final_table.txt")
dat <- dat[!(PH18)]

######################################################################
# Check link allele frequency and found in PH18
######################################################################
pdf("pdf/Figures/gDNA_allele_freq_scatterplots.pdf", 
    4, 
    4.5)
layout(matrix(1:2, ncol= 1), 
       heights = c(0.2, 1))
par(las= 2, 
    mgp= c(2, 0.5,0),
    tcl= -0.2)
dat[, {
  par(mar= c(0.2, 4.1, 1, 2.1))
  plot(density(NORMAL_alt_freq), 
       frame= F, 
       xaxt= "n", 
       main= paste(cdition, alt_class))
  par(mar= c(5.1, 4.1, 0.2, 2.1))
  plot(NORMAL_alt_freq, TUMOR_alt_freq, 
       col= adjustcolor(col, 0.8),
       pch= 16,
       cex= 0.8)
  leg <- .SD[, .N, .(class, col)]
  legend("topright",
         legend= leg[, paste0(class, " (" , N, ")")],
         fill= leg$col,
         bty= "n")
  abline(0, 1)
}, .(cdition, alt_class)]
dev.off()

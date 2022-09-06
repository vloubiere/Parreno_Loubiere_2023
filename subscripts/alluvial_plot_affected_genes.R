setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(cdition, system, FC_file)]
dat[, cdition:= factor(cdition, 
                       levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dat <- dcast(dat, 
             system+FBgn~cdition, 
             value.var = "diff")

pdf("pdf/Figures/alluvial_plot_timecourse_RNA.pdf", 
    width = 4, 
    height = 2.5)
par(cex= 0.7,
    las= 2,
    mar= c(4,5,1,6),
    tcl= -0.2,
    mgp= c(3,0.5,0),
    lwd= 0.5)
dat[, {
  pl <- vl_alluvial_plot(.SD[, -1],
                         col= c("cornflowerblue", "lightgrey", "tomato"),
                         ylab= "N genes")
  vars <- apply(.SD[, -1], 2, table)
  x <- rep(pl, each= 3)
  y <- c(apply(vars, 2, cumsum)-vars/2)
  text(x,
       y,
       c(vars),
       xpd= T,
       cex= 1)
}, system]
dev.off()

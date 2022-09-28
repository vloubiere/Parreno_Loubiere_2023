setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(FC_file, cdition, system)]
dat[, cdition:= factor(cdition, 
                       levels= c("PH18", "PHD11", "PHD9", "PH29"))]
setorderv(dat, 
          c("system", "cdition"))
dat <- dat[diff!="unaffected"]

pdf("pdf/RNA_boxplot_FC_up_down_genes.pdf", 2.8, 4)
par(mar= c(5,4,1,1),
    mgp= c(2,0.5,0),
    tcl= -0.2,
    las= 2)
dat[, {
  vl_boxplot(log2FoldChange~diff+cdition,
             col= c("cornflowerblue", "tomato"), 
             tilt.names= T, 
             ylab= "log2FoldChange",
             main= system)
  legend(c("topright", "bottomleft")[.GRP],
         c("down", "up"),
         fill= c("cornflowerblue", "tomato"),
         bty= "n")
  print("")
}, system]
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA"]
meta <- na.omit(meta[, .(cdition= gsub("^RNA_", "", cdition), FC_file)])
meta[, cdition:= factor(cdition, 
                        levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dat <- meta[, fread(FC_file), (meta)]
dat <- dcast(dat, 
             FBgn~cdition, 
             value.var = "diff")

pdf("pdf/Figures/alluvial_plot_timecourse_RNA.pdf", 
    width = 4, 
    height = 4)
par(cex= 0.7,
    las= 2,
    mar= c(2,2,2,1))
pl <- vl_alluvial_plot(dat[, -1],
                       class_levels = c("down", "unaffected", "up"),
                       col= c("cornflowerblue", "lightgrey", "tomato"),
                       ylab= "N genes")
vars <- apply(dat[, -1], 2, table)
x <- rep(pl, each= 3)
y <- c(apply(vars, 2, cumsum)-vars/2)
text(x,
     y,
     c(vars),
     xpd= T,
     cex= 0.6)
dev.off()

# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- fread("db/disc_sizes/discs_area_zfh1.txt")
setnames(dat, gsub("_", " ", gsub("__", "; ", names(dat))))
dat <- melt(dat, measure.vars = names(dat))
dat <- na.omit(dat)
dat[, counts:= paste0(" (n=" ,.N, ")"), variable]
dat[, variable:= paste0(variable, counts)]
dat[, variable:= factor(variable, unique(variable))]

pdf("pdf/Figure_2_ED_sizes_zfh1.pdf", 
    width = 2.25, 
    height = 3)
par(las= 2,
    tcl= -0.2,
    mgp= c(2.25,0.5,0),
    mar= c(6.5,4,1,0.5),
    lwd= 0.75)
vl_boxplot(value~variable,
           dat, 
           violin = T, 
           viocol = adjustcolor("rosybrown1", 0.6), 
           col= "white",
           compute_pval = list(c(1,2), c(2,4), c(1,3), c(3,4)), 
           viowex = 0.6,
           ylab= "Eye discs area (um2)", 
           tilt.names = T,
           ylim= c(0, 500000))
dev.off()
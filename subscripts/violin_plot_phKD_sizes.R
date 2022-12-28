# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- fread("db/disc_sizes/discs_area_ph_RNAi.txt")
dat <- dat[, grep("^PH", names(dat)), with= F]
setnames(dat, 
         c("PH18D", "PH29D", "PHD9", "PHD11"),
         c("no ph-KD", "Constant ph-KD ", "Transient ph-KD d9", "Transient ph-KD d11"))

pdf("pdf/Figure_1_ED_sizes.pdf", width = 2.25, height = 3.75)
par(las= 2,
    tcl= -0.2,
    mgp= c(2.25,0.5,0),
    mar= c(6.5,4,0.5,0.5),
    lwd= 0.75)
vl_boxplot(dat/1000, 
           violin = T, 
           viocol = adjustcolor("rosybrown1", 0.6), 
           col= "white",
           compute_pval = list(c(1,2), c(1,3), c(3,4), c(1,4)), 
           viowex = 0.6,
           ylab= "Eye discs area (mm2)", 
           tilt.names = T,
           ylim= c(0, 400))
dev.off()
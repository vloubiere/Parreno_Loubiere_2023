# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- fread("db/disc_sizes/discs_area_psc_RNAi.txt")
setnames(dat, 
         c("PSC18", "PSC29", "PSC 48H L1 KD"),
         c("No Psc/Su(z)2-KD", "Constant Psc/Su(z)2-KD ", "Transient Psc/Su(z)2-KD"))

pdf("pdf/Psc_RNAi_ED_sizes.pdf", width = 2.25, height = 3.75)
par(las= 2,
    tcl= -0.2,
    mgp= c(2.25,0.5,0),
    mar= c(6.5,4,4,0.5),
    lwd= 0.75)
vl_boxplot(dat, 
           violin = T, 
           viocol = adjustcolor("rosybrown1", 0.6), 
           col= "white",
           compute_pval = list(c(1,2), c(1,3), c(2,3)), 
           viowex = 0.6,
           ylab= "Eye discs area (um2)", 
           tilt.names = T,
           srt= 30)
dev.off()
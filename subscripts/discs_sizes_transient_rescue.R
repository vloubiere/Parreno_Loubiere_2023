setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- fread("db/disc_sizes/discs_area_rescue_transient_RNAi.txt")
dat <- dat[, .(`gfpRNAi, wRNAi`, `stat92eRNAi, wRNAi`, `zfh1RNAi, wRNAi`, `gfpRNAi, phRNAi`, `stat92eRNAi, phRNAi`, `zfh1RNAi, phRNAi`)]
setnames(dat,
         gsub("RNAi, ", " + ", names(dat)))
setnames(dat,
         gsub("stat92e", "Stat92E", names(dat)))
setnames(dat,
         gsub("RNAi", "-KD", names(dat)))
setnames(dat, 
         paste0(names(dat), " (n=", apply(dat, 2, function(x) sum(!is.na(x))), ")"))

# Plot ----
pdf("pdf/review_disc_sizes_transient_RNAi.pdf", 3, 3)
vl_par(mgp= c(1.7, 0.2, 0))
vl_boxplot(dat, 
           violin = T, 
           viocol = adjustcolor("lightgrey", 0.6), 
           col= "white",
           compute.pval = list(c(1,2), c(1,3), c(1,4), c(4,5), c(4,6)), 
           viowex = 0.6,
           ylab= "Eye disc area (um2)", 
           tilt.names = T,
           ylim= c(0, 500000),
           lwd= .75)
dev.off()
setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- fread("db/disc_sizes/discs_area_rescue_constant_RNAi.txt",
             sel= c("gfpRNAi, wRNAi",
                    "stat92eRNAi, wRNAi",
                    "zfh1RNAi, wRNAi",
                    "gfpRNAi, phRNAi",
                    "stat92eRNAi, phRNAi",
                    "zfh1RNAi, phRNAi"))
setnames(dat,
         function(x) gsub("stat92e", "Stat92E", x))
setnames(dat,
         function(x) gsub("RNAi, ", "+", x))
setnames(dat,
         function(x) gsub("RNAi$", " RNAi", x))
dat <- melt(dat, measure.vars = names(dat))
dat <- na.omit(dat)
dat[, variable:= paste0(variable, " (", .N, ")"), variable]
dat[, variable:= factor(variable, unique(variable))]

pdf("pdf/review_rescue_constant_RNAi.pdf", 3.25, 2.75)
vl_par(mgp= c(2,0.35,0),
       lwd= 0.75)
vl_boxplot(value~variable,
           dat, 
           violin = T, 
           viocol = adjustcolor("lightgrey", 0.6), 
           col= "white",
           compute.pval = list(c(1,2), c(1,3), c(1,4), c(4,5), c(4,6)), 
           viowex = 0.6,
           ylab= "Eye discs area (um2)", 
           tilt.names = T,
           ylim= c(0, 500000))
dev.off()
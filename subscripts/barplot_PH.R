require(vlfunctions)

# Import
dat <- fread("Rdata/final_gene_features_table.txt")
# Extend promoters
dat <- vl_resizeBed(dat, "start", 250, 250)

# Overlaps
PH <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
dat[, PRC1_PH18:= vl_covBed(.SD, PH)>0]
PH <- vl_importBed("db/peaks/cutnrun/PH_PH29_confident_peaks.narrowPeak")
dat[, PRC1_PH29:= vl_covBed(.SD, PH)>0]
PH <- vl_importBed("db/peaks/cutnrun/PH_PHD11_confident_peaks.narrowPeak")
dat[, PRC1_PHD11:= vl_covBed(.SD, PH)>0]

# Prepare for ploting
pl <- melt(dat, "recovery", patterns("^PRC1_PH"))
pl[, recovery:= switch(recovery, 
                       "noRecovery"= "Irreversible",
                       "Recovery"= "Reversible"), recovery]
pl[, bound:= sum(value), .(variable, recovery)]
pl[, unbound:= .N-sum(value), .(variable, recovery)]
pl <- unique(pl[!is.na(recovery), .(recovery, variable, bound, unbound)])
setorderv(pl, c("recovery", "variable"))
pl[, cdition:= paste0(recovery, "_", variable)]
mat <- as.matrix(pl[, .(cdition, bound, unbound)], 1)

pdf("pdf/Figure_3_barplots_PH_binding.pdf", width= 5, height = 3.5)
par(las= 1,
    tcl= -0.2,
    mgp= c(2, 0.5, 0),
    mar= c(5,5,2,10))
bar1 <- barplot(t(mat)[,1:3],
               width = 1,
               space = c(0, 0.25, 0.25),
               xaxt= "n",
               ylab= "Number of gene promoters",
               col= c("rosybrown1", "white"),
               xlim= c(0,7))
barplot(t(mat)[,1:3],
        width = 1,
        space = c(0, 0.25, 0.25),
        xaxt= "n",
        col= c("grey40", 
               adjustcolor("rosybrown1", 0.6)),
        xlim= c(0,7), 
        density = c(40, -1),
        add= T)
bar2 <- barplot(t(mat)[,4:6],
               width = 1,
               space = c(4,0.25,0.25),
               xaxt= "n",
               col= c("palegreen3","white"),
               xlim= c(0,7),
               add= T)
bar2 <- barplot(t(mat)[,4:6],
                width = 1,
                space = c(4,0.25,0.25),
                xaxt= "n",
                col= c("grey40",
                       adjustcolor("palegreen3", 0.6)),
                xlim= c(0,7),
                density = c(40, -1),
                add= T)
vl_tilt_xaxis(c(bar1, bar2), 
              labels= rep(c("No ph-KD", "Constant ph-KD", "Transient ph-KD"), 2), 
              srt= 30)
legend(par("usr")[2],
       par("usr")[4],
       fill= c(adjustcolor("rosybrown1", 0.6), 
               adjustcolor("palegreen3", 0.6),
               "black"),
       density= c(-1,-1,40),
       legend= c("Irreversible", "Reversible", "PH-bound promoter"),
       bty= "n",
       xpd= NA)
dev.off()
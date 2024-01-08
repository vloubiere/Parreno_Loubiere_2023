setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- readRDS("Rdata/final_gene_features_table.rds")[class %in% c("Irreversible", "Reversible")]
dat[, class:= droplevels(class)]
dat <- dat[, .(class, seqnames, start, end, strand, FBgn, symbol)]

# K27Ac Binding ----
id <- vl_closestBed(dat, "db/peaks/cutnrun/H3K27Ac_PH18_confident_peaks.narrowPeak")[between(dist, -2500, 0), FBgn]
dat[, `Ac_18`:= FBgn %in% id]
id <- vl_closestBed(dat, "db/peaks/cutnrun/H3K27Ac_PH29_confident_peaks.narrowPeak")[between(dist, -2500, 0), FBgn]
dat[, `Ac_29`:= FBgn %in% id]
id <- vl_closestBed(dat, "db/peaks/cutnrun/H3K27Ac_PHD11_confident_peaks.narrowPeak")[between(dist, -2500, 0), FBgn]
dat[, `Ac_D11`:= FBgn %in% id]

.m <- melt(dat,
           id.vars = "class",
           measure.vars = c("Ac_18", "Ac_29", "Ac_D11"))
pl <- .m[, .(value= sum(value), total= .N), .(variable, class)]

# Plot ----
Cc <- rep(c("rosybrown1", "palegreen3"), each= 3)
Cc <- adjustcolor(Cc, .6)

pdf("pdf/review_K27Ac_binding_irreversible_reversible.pdf", 3, 3)
vl_par(lend= 2,
       mgp= c(1,.25, 0),
       lwd= 0.75)
bar <- barplot(total~variable+class,
               pl, 
               beside= T,
               col= Cc,
               ylab= "Number of genes",
               main= "K27Ac binding",
               xaxt= "n",
               xlab= NA,
               width= 1,
               space= c(0, 0.2))
vl_tilt_xaxis(bar,
              labels = rep(c("No ph-KD", "Constant ph-KD", "Transient ph-KD"), 2),
              srt= 25)
par(lwd= 0.5)
barplot(value~variable+class,
        pl,
        beside= T,
        add= T,
        dens= 40,
        col= "black",
        axes= F,
        xaxt= "n",
        yaxt= "n",
        width= 1,
        space= c(0, 0.2))
legend(par("usr")[2],
       par("usr")[4],
       fill= c(unique(Cc), "black"),
       legend= c("Irreversible", "Reversible", "Bound"),
       dens= c(NA, NA, 40),
       bty= "n",
       xpd= NA,
       cex= .6)
dev.off()
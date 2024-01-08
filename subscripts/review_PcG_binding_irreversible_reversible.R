setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[class %in% c("Irreversible", "Reversible")]
dat[, class:= droplevels(class)]
dat <- dat[, .(class, seqnames, start, end, strand, FBgn, symbol, PcG_bound)]
setnames(dat, "PcG_bound", "PcG_bound_18")

# PH29 ----
ov <- vl_closestBed(dat, "db/bed/merged_K27_domains/PH29.bed")
ov <- ov[dist==0]
ov[, max_start:= apply(.SD, 1, max), .SDcols= c("start", "start.b")]
ov[, min_end:= apply(.SD, 1, max), .SDcols= c("end", "end.b")]
ov <- ov[(min_end-max_start+1)/(end-start+1)>0.5] # 50% of the gene body covered with the mark
dat[, PcG_bound_29:= FBgn %in% ov$FBgn]

# PHD11 ----
ov <- vl_closestBed(dat, "db/bed/merged_K27_domains/PHD11.bed")
ov <- ov[dist==0]
ov[, max_start:= apply(.SD, 1, max), .SDcols= c("start", "start.b")]
ov[, min_end:= apply(.SD, 1, max), .SDcols= c("end", "end.b")]
ov <- ov[(min_end-max_start+1)/(end-start+1)>0.5] # 50% of the gene body covered with the mark
dat[, PcG_bound_D11:= FBgn %in% ov$FBgn]

# melt
.m <- melt(dat,
           id.vars = "class",
           measure.vars = c("PcG_bound_18", "PcG_bound_29", "PcG_bound_D11"))
pl <- .m[, .(value= sum(value), total= .N), .(variable, class)]

# Plot ----
Cc <- rep(c("rosybrown1", "palegreen3"), each= 3)
Cc <- adjustcolor(Cc, .6)

pdf("pdf/review_PcG_binding_irreversible_reversible.pdf", 3, 3)
vl_par(lend= 2,
       mgp= c(1,.25, 0),
       lwd= 0.75)
bar <- barplot(total~variable+class,
               pl, 
               beside= T,
               col= Cc,
               ylab= "Number of genes",
               main= "PcG binding",
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
       legend= c("Irreversible", "Reversible", "Within K27me3 domain"),
       dens= c(NA, NA, 40),
       bty= "n",
       xpd= NA,
       cex= .6)
dev.off()
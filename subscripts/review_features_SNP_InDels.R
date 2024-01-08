setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- readRDS("Rdata/review_gDNA_final_table.rds")
dat <- dat[allele_freq>0.2]
dat[, class:= fcase(any(grepl("^no ph-KD", cdition)), "anyCtl",
                    tumor_fraction %in% levels(tumor_fraction)[-1], "atLeastTwoTum"), .(id, tumor_fraction)]
dat <- unique(dat[!is.na(class), .(id, class, annotation)])

# prepare plotting
feat <- table(dat$annotation, dat$class)
perc <- apply(feat, 2, function(x) x/sum(x)*100)
Cc <- rainbow(length(unique(dat$annotation)))
Cc <- adjustcolor(Cc, 0.3)

# Plot
pdf("pdf/review_SNP_InDel_features.pdf", 3, 3)
vl_par(mai= c(.9, .9, .9, 1.2))
bar <- barplot(perc, 
               col= Cc,
               ylab= "SNV/InDels (%)",
               border= NA,
               xaxt= "n")
vl_tilt_xaxis(bar,
              labels = colnames(perc))
text(rep(bar, each= nrow(feat)),
     apply(perc, 2, cumsum)-perc/2,
     ifelse(perc>7, feat, NA),
     cex= 0.5)
legend(par("usr")[2],
       par("usr")[4],
       fill= Cc,
       legend= rownames(feat),
       xpd= T,
       bty= "n",
       cex= 0.5, 
       y.intersp = 0.8,
       border= NA)
dev.off()
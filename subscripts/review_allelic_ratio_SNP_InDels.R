setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- readRDS("Rdata/review_gDNA_final_table.rds")

# Plot ----
pdf("pdf/review_SNP_InDel_allelic_fraction.pdf", 3, 3)
vl_par(mai= c(1.2, .9, 1.2, .9))
dens <- density(dat$allele_freq)
plot(dens,
     type= "n",
     frame= F,
     xlab= "SNV/InDel allele frequency",
     main= NA,
     xaxt= "n")
axis(1, padj = -1.25, gap.axis = 0)
axis(2, gap.axis = 0)
polygon(dens$x,
        dens$y,
        border= NA,
        col ="lightgrey")
abline(v= 0.2,
       lty= "11")
tab <- table(dat$allele_freq>0.2)
text(c(0.1, 0.6),
     par("usr")[4],
     paste0(round(tab/sum(tab)*100, 1), "%\n", formatC(tab, big.mark = ",")),
     xpd= T,
     pos= 3,
     cex= .3)
dev.off()
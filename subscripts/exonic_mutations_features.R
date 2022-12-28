# setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!is.na(exonic_annotation)]
dat[, exonic_annotation:= factor(exonic_annotation)]
dat[, exonic_annotation:= droplevels(exonic_annotation)]
dat[, cdition:= switch(cdition,
                       "PH18_2"= "no ph-KD",
                       "PH29_1"= "Constant ph-KD.1",
                       "PH29_2"= "Constant ph-KD.2",
                       "PHD9_1"= "Trans. ph-KD d9.1",
                       "PHD9_2"= "Trans. ph-KD d9.2",
                       "PHD11_1"= "Trans. ph-KD d11.1",
                       "PHD11_2"= "Trans. ph-KD d11.2"), cdition]
dat[, cdition:= factor(cdition, 
                       c("Trans. ph-KD d11.2",
                         "Trans. ph-KD d11.1",
                         "Trans. ph-KD d9.2",
                         "Trans. ph-KD d9.1",
                         "Constant ph-KD.2",
                         "Constant ph-KD.1",
                         "no ph-KD"))]
Cc <- rainbow(length(levels(dat$exonic_annotation)))
Cc <- adjustcolor(Cc, 0.3)

# Format before plotting
pl <- dcast(dat, 
            exonic_annotation~cdition, 
            value.var = "id", 
            fun.aggregate = length, 
            fill = 0)
mat <- as.matrix(pl, 1)
perc <- apply(mat, 2, function(x) x/sum(x)*100)


pdf("pdf/Figure_1_exonic_mutations_features.pdf", 5, 2.75)
par(mar= c(4,8,1,9),
    las= 1,
    lwd= 0.1,
    tcl= -0.2,
    mgp= c(1,0.15, 0),
    cex.axis= 9/12)
bar <- barplot(perc, 
               col= Cc,
               xlab= "Exonic mutations (%)",
               border= NA,
               horiz= T)
mat[perc<10] <- NA
text(apply(perc, 2, cumsum)-perc/2,
     rep(bar, each= nrow(mat)),
     mat,
     cex= 0.5)
legend(par("usr")[2],
       par("usr")[4],
       fill= Cc,
       legend= levels(dat$exonic_annotation),
       xpd= T,
       bty= "n",
       cex= 0.5, 
       y.intersp = 0.8,
       border= NA)
dev.off()
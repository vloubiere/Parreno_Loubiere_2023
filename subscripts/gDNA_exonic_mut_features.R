setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[(cdition=="PH18_2" | !(PH18)) & type!="unknown"]
dat[cdition=="PH18_2", occurence:= "Ph18_control"]
dat <- unique(dat[, .(id, type, class, occurence)])
dat[, type:= factor(type)]
dat[, type:= droplevels(type)]
dat[, occurence:= factor(occurence,
                         c("Ph18_control", 
                           "single condition",
                           "shared >=1 conditions"))]
Cc <- rainbow(length(levels(dat$type)))
Cc <- adjustcolor(Cc, 0.4)

pdf("pdf/gDNA_exonic_mut_features.pdf", 
    width = 4, 
    height = 4.5)
par(mfrow= c(1,3),
    mar= c(15,3,4,2),
    las= 2,
    lwd= 0.1,
    tcl= -0.2,
    mgp= c(2,0.5, 0))
dat[, {
  pl <- .SD[, .N, .(type, occurence)]
  pl <- dcast(pl, type~occurence, value.var = "N", fill = 0)
  mat <- as.matrix(pl, 1)
  perc <- apply(mat, 2, function(x) x/sum(x)*100)
  bar <- barplot(perc, 
                 main= class[1], 
                 col= Cc,
                 ylab= "Percentage (%)")
  mat[perc<2] <- NA
  text(rep(bar, each= nrow(mat)),
       apply(perc, 2, cumsum)-perc/2,
       mat,
       cex= 0.5)
  if(.GRP%%2==0)
  {
    plot.new()
    legend(grconvertX(0, "nfc", "user"),
           par("usr")[4],
           fill= Cc,
           legend= levels(type),
           xpd= T,
           bty= "n",
           cex= 0.8)
  }
  print(class)
}, class]
dev.off()


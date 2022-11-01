setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[cdition=="PH18_2" | !(PH18)]
dat[cdition=="PH18_2", occurence:= "PH18_control"]
dat <- unique(dat[, .(id, annotation, class, occurence)])
dat[, annotation:= factor(annotation)]
dat[, annotation:= droplevels(annotation)]
dat[occurence=="single condition", occurence:=  "single condition (not PH18)"]
dat[occurence=="shared >=1 conditions", occurence:=  "shared >=1 conditions (not PH18)"]
dat[, occurence:= factor(occurence,
                         c("PH18_control", 
                           "single condition (not PH18)",
                           "shared >=1 conditions (not PH18)"))]
Cc <- rainbow(length(levels(dat$annotation)))
Cc <- adjustcolor(Cc, 0.4)

pdf("pdf/gDNA_mut_features.pdf", 
    width = 4, 
    height = 4.5)
par(mfrow= c(1,3),
    mar= c(15,3,4,2),
    las= 2,
    lwd= 0.1,
    tcl= -0.2,
    mgp= c(2,0.5, 0))
dat[, {
  pl <- dcast(.SD, 
              annotation~occurence, 
              value.var = "id", 
              fun.aggregate = length, 
              fill = 0)
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
           legend= levels(annotation),
           xpd= T,
           bty= "n",
           cex= 0.8)
  }
  print(class)
}, class]
dev.off()


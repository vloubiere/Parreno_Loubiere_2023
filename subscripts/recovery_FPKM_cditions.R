setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)
require(preprocessCore)
require(basicPlotteR)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[, col:= fcase(recovery=="Recovery", "palegreen3", 
                  recovery=="noRecovery", "rosybrown1",
                  default= "grey")]
cols <- names(dat[, FPKM_PH18:FPKM_PHD9])
# quantile normalize FPKMs
dat[, mark:= symbol %in% c("Ets21C", 
                           "Sox15", 
                           "lncRNA:iab8", 
                           "upd1", 
                           "upd2", 
                           "upd3", 
                           "en", 
                           "wg", 
                           "eve", 
                           "inv", 
                           "fd102c")]
.m <- melt(dat,
           id.vars = c("symbol", "col", "recovery", "mark"),
           measure.vars = c("FPKM_PH18", "FPKM_PH29", "FPKM_PHD11", "FPKM_PHD9"), 
           value.name = "FPKM")
.m[, percentile:= frank(FPKM)/.N*100, variable]
Cc <- c("rosybrown1", "palegreen3")
Cc <- adjustcolor(Cc, 0.7)

pdf("pdf/recovery_groups_FPKMs.pdf", 8.5, 2)
par(mar= c(3.5,4,2,1),
    mgp= c(1.75,0.5,0),
    tcl= -0.2,
    las= 2,
    mfrow= c(1,4))
# Boxplot FPKM
vl_boxplot(FPKM~variable+recovery, 
           .m, 
           las=2,
           col= rep(Cc, each= 4),
           compute_pval= list(c(1,2), c(1,3), c(1,4), c(5,6), c(5,7), c(5,8)),
           names= function(x) gsub(".Recovery|.noRecovery|FPKM_", "", x),
           tilt.names= T,
           ylab= "FPKM")
legend("topleft",
       bty= "n",
       legend= c("noRecovery", "Recovery"),
       fill= Cc,
       inset = c(0, -0.2),
       xpd= NA)
# Boxplot percentile
vl_boxplot(percentile~variable+recovery, 
           .m, 
           las=2,
           col= rep(Cc, each= 4),
           compute_pval= list(c(1,2), c(1,3), c(1,4), c(5,6), c(5,7), c(5,8)),
           names= function(x) gsub(".Recovery|.noRecovery|FPKM_", "", x),
           tilt.names= T,
           ylab= "FPKM percentile")
# FPKM examples
plot(NA,
     xlim= c(0.6,5.5),
     ylim= c(-10,7.5),
     ylab= "log2(FPKM+0.001)",
     xaxt= "n",
     xlab= NA)
axis(1, 
     at= 1:4, 
     labels= c("PH18", "PH29", "PHD11", "PHD9"))
dat[(mark), {
  x <- 1:4
  y <- log2(unlist(.SD)+0.001)
  col <- Cc[ifelse(recovery=="noRecovery", 1, 2)]
  lines(x, y, 
        lwd= 2,
        col= col)
  points(x, y, 
        pch= 16,
        col= col,
        cex= 1.5)
  text(4, y[4], symbol[1], col= col, pos= 4)
}, FBgn, .SDcols= cols]
# Scaled FPKMs examples
plot(NA,
     xlim= c(0.6,5.5),
     ylim= c(0,16),
     ylab= "Scaled log2(FPKM+0.001)",
     xaxt= "n",
     xlab= NA)
axis(1, 
     at= 1:4, 
     labels= c("PH18", "PH29", "PHD11", "PHD9"))
dat[(mark), {
  x <- 1:4
  y <- log2(unlist(.SD)+0.001)
  y <- y-min(y)
  col <- Cc[ifelse(recovery=="noRecovery", 1, 2)]
  lines(x, 
        y, 
        lwd= 2,
        col= col)
  points(x, 
         y, 
         pch= 16,
         col= col,
         cex= 1.5)
  text(4, y[4], symbol[1], col= col, pos= 4)
}, FBgn, .SDcols= cols]
dev.off()

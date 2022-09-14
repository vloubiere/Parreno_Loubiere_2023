setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[, col:= fcase(recovery=="Recovery", "palegreen3", 
                  recovery=="noRecovery", "rosybrown1",
                  default= "grey")]

.m <- melt(dat,
           id.vars = c("symbol", "col", "recovery"),
           measure.vars = c("FPKM_PH18", "FPKM_PH29", "FPKM_PHD11", "FPKM_PHD9"), 
           value.name = "FPKM")
.m[, FPKM:= log2(FPKM+0.001)]
setorderv(.m, "FPKM")
.m[, x:= seq(.N), variable]
.m[symbol %in% c("Ets21C", "Sox15", "lncRNA:iab8", "upd1", "upd2", "upd3", 
                 "en", "wg", "eve", "inv", "fd102c"), mark:= T]
.m[, percentile:= x/.N*100, variable]
leg <- na.omit(unique(dat[, .(recovery, col)]))

pdf("pdf/Figures/recovery_groups_FPKMs.pdf", 8.5, 2)
par(mar= c(2,4,2,1),
    mgp= c(1.75,0.5,0),
    tcl= -0.2,
    las= 1,
    mfrow= c(1,4))
.m[, {
  plot(x,
       FPKM,
       type= "l",
       lwd= 2,
       pch= 16,
       col= "grey",
       main= variable,
       xaxt= "n",
       xlab= NA,
       ylim= c(-11, 13.5))
  title(xlab= "Rank", line= 0.2)
  leg[, legend("topleft",
               fill= col,
               legend= recovery,
               bty= "n",
               border= NA,
               cex= 0.8)]
  .SD[!is.na(recovery), {
    points(x,
           FPKM,
           pch= 16,
           col= col) 
  }]
  .SD[(mark), {
    points(x,
           FPKM,
           pch= 16,
           col= col)
    segments(x,
             FPKM,
             x+strwidth("M"),
             FPKM-strheight("M"),
             col= col)
    text(x+strwidth("M")*1.1,
         FPKM-strheight("M")*1.5,
         symbol,
         pos= 4,
         offset= 0,
         col= col)
    
  }]
}, variable]
par(mar= c(2,4,1.5,1))
vl_boxplot(percentile~variable+recovery, 
           .m, 
           names= function(x) gsub(".Recovery|.noRecovery|FPKM_", "", x),
           tilt.names= T,
           ylab= "FPKM percentile",
           col= rep(c("rosybrown1", "palegreen3"), each= 4))
dev.off()







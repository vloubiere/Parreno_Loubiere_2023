setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[!is.na(class)]
setorderv(dat, "class")
dat[, class:= paste0(class, " (", formatC(.N, big.mark = ","), ")"), class]
dat[, class:= factor(class, unique(class))]

# Format ----
FC <- melt(dat,
           id.vars = "class",
           measure.vars = c("log2FoldChange_PH18", "log2FoldChange_PH29", "log2FoldChange_PHD11"))

FPKM <- melt(dat,
             id.vars = "class",
             measure.vars = c("FPKM_PH18", "FPKM_PH29", "FPKM_PHD11"))

# Ploting function ----
.f <- function(x, ylab)
{
  Cc <- c("lightgrey", "rosybrown1", "palegreen3")
  Cc <- adjustcolor(Cc, .6)
  vl_boxplot(value~class+variable,
             data= x,
             ylab= ylab,
             tilt.names= T, 
             col= Cc, 
             boxwex= 0.6,
             xaxt= "n",
             at= c(1,2,3,6,7,8,11,12,13),
             compute.pval= list(c(1,2), c(2,3), c(1,3),
                                c(4,5), c(5,6), c(4,6),
                                c(7,8), c(8,9), c(7,9)),
             pval.cex= 5/12,
             lwd= 0.5)
axis(1,
     at= c(2, 7, 12),
     labels = NA)
cditions <- c("No ph-KD", "Constant ph-KD", "Transient ph-KD")
vl_tilt_xaxis(x= c(2, 7, 12),
              y = grconvertY(grconvertY(0, "npc", "inch")-strheight("M", "inch"), "inch", "user"),
              labels = cditions,
              srt= 25)
legend("topright",
       fill= Cc,
       legend= levels(x$class),
       bty= "n",
       xpd= T,
       cex= 7/12,
       inset= c(-0.6, -0.2))
}

# Plot ----
pdf("pdf/review_RNAseqFC_FPKM_per_class.pdf",
    width = 3, 
    height = 3)
vl_par(lwd= 0.5)
.f(FC, "RNA-Seq fold change (log2)")
.f(FPKM, "RNA-Seq FPKM")
dev.off()
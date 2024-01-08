setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- fread("db/DNA_repair/pH2Av_foci_per_cell.txt")
pl <- melt(dat, measure.vars = names(dat))
pl <- na.omit(pl)
pl <- pl[, .(values= .(value),
             mean= mean(value),
             se= sqrt(sum((value-mean(value))^2/(length(value)-1)))/sqrt(length(value))), variable]


# Plot ----
Cc <- c("grey80", "grey40")

pdf("pdf/review_DNA_repair_barplot.pdf", 3, 3)
vl_par(xpd= NA,
       mgp= c(1, .2,0))
pl[, {
  bar <- barplot(mean,
                 col= rep(Cc, each= 3),
                 space = c(1,0.2,0.2,1,0.2,0.2),
                 ylab= "pH2Av foci/cell")
  axis(1, bar, rep(c("0'", "30'", "480'"), 2), lwd= 0, gap.axis= 0, padj= -1.25)
  arrows(bar, mean, bar, mean-se, angle = 90, length = 0.025)  
  arrows(bar, mean, bar, mean+se, angle = 90, length = 0.025)
  
  y <- par("usr")[3]-strheight("M")*1.25
  for(cmb in list(c(1, 3, "No ph-KD"), c(4, 6, "Trans. ph-KD")))
  {
    x0 <- bar[as.numeric(cmb[1])]
    x1 <- bar[as.numeric(cmb[2])]
    segments(x0, y, x1, y)
    text(mean(c(x0, x1)),
         y,
         cmb[3],
         cex= 8/12,
         pos= 1,
         offset= 0.25)
  }
  
  # Test control ----
  for(cmb in list(c(1,4), c(2,3), c(5,6)))
  {
    y <- max(mean[cmb]+se[cmb])+strheight("M")/3
    if(identical(cmb, c(1,4)))
      y <- y+12.15
    x0 <- bar[cmb[1]]
    x1 <- bar[cmb[2]]
    segments(x0, y, x1, y)
    pval <- t.test(values[[cmb[1]]], values[[cmb[2]]])$p.value
    vl_plot_pval_text(mean(c(x0, x1)),
                      y+strheight("M", cex= 0.7),
                      pval,
                      stars.only = T,
                      cex= 7/12)
    FC <- mean[cmb[2]]/mean[cmb[1]]
    text(mean(c(x0, x1)),
         y,
         paste0("x", round(FC, 2)),
         pos= 3,
         offset= 0.1,
         cex= 5/12)
  }
}]
dev.off()
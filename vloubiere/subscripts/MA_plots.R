setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

dat <- readRDS("Rdata/RNA_tables_object.rds")$FC

lapply(seq(dat), function(i){
  lim <- 7
  .c <- dat[[i]]
  png(paste0("pdf/MA_plots/MA_Plots_", names(dat)[i], ".png"), 
      width = 1000, 
      height = 750)
  par(mfrow= c(3,4),
      pty= "s")
  .c[, 
     {
       y <- log2FoldChange
       pch <- ifelse(abs(y)<=lim, 16, 2)
       y[y >  lim] <- lim
       y[y < -lim] <- -lim
       Cc <- ifelse(is.na(padj) | padj >= 0.01, "lightgrey", ifelse(y>0, "red", "blue"))
       
       plot(baseMean, 
            y, 
            col= Cc, 
            ylim= c(-lim, lim), 
            pch= pch, 
            cex= 0.5, 
            log= "x", 
            ylab= "log2FoldChange", 
            xlab= "baseMean",
            las= 1)
       mtext(cdition, 
             line = 1, 
             cex= 0.7)
       leg <- c(paste(length(which(Cc=="red")), "Up"), 
                paste(length(which(Cc=="blue")), "Down"), "(padj<0.01)")
       legend("topright", 
              leg, 
              bty= "n", 
              text.col = c("red", "blue", "black"), 
              cex = 0.8)
       print("DONE")
     }, cdition]
  dev.off()
})


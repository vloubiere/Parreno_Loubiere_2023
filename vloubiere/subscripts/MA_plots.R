dat <- data.table(file= list.files("db/FC_tables/", "RNA_epiCancer", full.names = T))
dat[, cdition:= gsub("RNA_epiCancer_|.txt$" , "", basename(file)), file]
dat <- dat[, fread(file), (dat)]

dir.create("pdf/MA_plots", showWarnings = F)

dat[, 
    {
      y <- log2FoldChange
      pch <- ifelse(abs(y)<=5, 16, 2)
      y[y >  5] <- 5
      y[y < -5] <- -5
      Cc <- ifelse(is.na(padj) | padj >= 0.01, "lightgrey", ifelse(y>0, "red", "blue"))
      
      pdf(paste0("pdf/MA_plots/", cdition, ".pdf"), width = 7, height = 8)
      plot(baseMean, 
           y, 
           col= Cc, 
           ylim= c(-5, 5), 
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
      dev.off()

      print("DONE")
      
}, cdition]

setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)
require(pheatmap)

# Import
dat <- data.table(files= sapply(c("EzJ9_vs_WKD", "EzJ11_vs_WKD", "Ez29_vs_W29", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29"), function(x) 
  list.files("db/FC_tables/", x, full.names = T)))
dat[, cdition:= tstrsplit(files, "epiCancer_|_FC", keep= 2)]
dat <- dat[, fread(files, select = c(1,2,3,6), col.names = c("FBgn", "baseMean", "log2FoldChange", "padj")),  (dat)]

pdf("pdf/MA_plots_PH_vs_W.pdf", width = 7, height = 8)
par(mfrow= c(2, 2), las= 1)
dat[, 
    {
      y <- log2FoldChange
      pch <- ifelse(abs(y)<=5, 16, 2)
      y[y >  5] <- 5
      y[y < -5] <- -5
      Cc <- ifelse(is.na(padj) | padj >= 0.01, "lightgrey", ifelse(y>0, "red", "blue"))
      plot(baseMean, y, col= Cc, ylim= c(-5, 5), pch= pch, cex= 0.5, log= "x", ylab= "log2FoldChange", xlab= "baseMean")
      mtext(cdition, line = 1, cex= 0.7)
      leg <- c(paste(length(which(Cc=="red")), "Up"), paste(length(which(Cc=="blue")), "Down"), "(padj<0.01)")
      legend("topright", leg, bty= "n", text.col = c("red", "blue", "black"), cex = 0.8)
      print("DONE")
    }, cdition]
dev.off()

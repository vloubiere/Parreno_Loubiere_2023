setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[project=="RNA_epiCancer"]
meta <- meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), .(DESeq2_object, cdition)]
meta <- unique(meta[FC_file!="NA"])
meta <- meta[, fread(FC_file), (meta)]
meta[, cdition:= factor(cdition, 
                        levels= c("RNA_PH18", "RNA_PHD11", "RNA_PHD9", "RNA_PH29", "RNA_PHD11_T2", "RNA_PHD11_T3"))]
setorderv(meta, c("DESeq2_object", "cdition"))

lim <- 10

pdf("pdf/RNA/MA_plots.pdf", 10, 8.5)
par(mfrow=c(3,4))
meta[, {
  y <- log2FoldChange
  pch <- ifelse(abs(y)<=lim, 16, 2)
  y[y >  lim] <- lim
  y[y < -lim] <- -lim
  Cc <- ifelse(padj<=0.05 & abs(y)>1, ifelse(y>0, "red", "blue"), "lightgrey")
  
  pl_ord <- order(Cc=="lightgrey", # Color dots plot last 
                  decreasing = T)
  
  plot(baseMean[pl_ord], 
       y[pl_ord], 
       col= Cc[pl_ord], 
       ylim= c(-lim, lim), 
       pch= pch, 
       cex= 0.5, 
       log= "x", 
       ylab= "log2FoldChange", 
       xlab= "baseMean",
       las= 1)
  mtext(paste0(DESeq2_object, "\n", 
               gsub(paste0(DESeq2_object, "_(.*).txt$"), "\\1", basename(FC_file))), 
        line = 1, 
        cex= 0.7)
  leg <- c(paste(length(which(Cc=="red")), "Up"), 
           paste(length(which(Cc=="blue")), "Down"), "(p<0.05, |log2FC|>1)")
  legend("topright", 
         leg, 
         bty= "n", 
         text.col = c("red", "blue", "black"), 
         cex = 0.5)
  print("DONE")
}, .(DESeq2_object, cdition, FC_file)]
dev.off()
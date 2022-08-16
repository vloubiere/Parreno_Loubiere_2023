setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[project=="RNA_epiCancer"]
meta <- unique(meta[FC_file!="NA", .(DESeq2_object, cdition, FC_file)])
dat <- meta[, fread(FC_file), (meta)]
dat[, cdition:= factor(cdition, 
                        levels= c("RNA_PH18", "RNA_PHD11", "RNA_PHD9", "RNA_PH29"))]
setorderv(dat, c("DESeq2_object", "cdition"))
dat[, col:= switch(diff, "up"= "red", "down"= "blue", "unaffected"= "lightgrey"), diff]
dat <- dat[order(col=="lightgrey", decreasing= T)]

pdf("pdf/Figures/MA_plots_RNA.pdf", 10, 8.5)
par(mfrow=c(3,4))
dat[, {
  plot(baseMean, 
       ifelse(abs(log2FoldChange)>10, sign(log2FoldChange)*10, log2FoldChange),
       col= col,
       ylim= c(-10, 10),
       pch= ifelse(abs(log2FoldChange)>10, 17, 16),
       cex= 0.5,
       log= "x", 
       ylab= "log2FoldChange", 
       xlab= "baseMean",
       las= 1)
  title(main= paste(cdition,
                    fcase(grepl("GFP+", fixed = T, DESeq2_object), "GFP+",
                          grepl("GFP-", fixed = T, DESeq2_object), "GFP-")))
  leg <- c(paste(sum(col=="red"), "Up"), 
           paste(sum(col=="blue"), "Down"), "(p<0.05, |log2FC|>1)")
  legend("topright", 
         leg, 
         bty= "n", 
         text.col = c("red", "blue", "black"), 
         cex = 0.5)
  print("DONE")
}, .(DESeq2_object, cdition, FC_file)]
dev.off()

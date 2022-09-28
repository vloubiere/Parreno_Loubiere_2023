setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
bins <- vl_binBSgenome("dm6", 
                       bins_width = 5e3, 
                       steps_width = 5e3, 
                       restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"))

# cutnrun
dat <- meta[!is.na(bw_file), .(score= vl_bw_coverage(bins, bw_file)), .(ChIP, cdition, rep)]
dat[, idx:= seq(.N), .(ChIP, cdition, rep)]
dat <- dcast(dat,
             idx~ChIP+cdition+rep,
             value.var = "score")
mat <- as.matrix(dat, 1)

# Plot
pdf("pdf/cutnrun_PCC_replicates.pdf",
    width = 12,
    height = 5.5)
par(las= 1,
    tcl= -0.2,
    mgp= c(2,0.5,0))
for(patt in c("^H", "^PH_"))
{
  if(patt=="^H")
    par(mar= c(5.1, 4.1, 4.1, 2.1)) else
      par(mar= c(5.1, 25, 4.1, 25)) 
  sub <- mat[, grep(patt, colnames(mat))]
  cditions <- gsub("_rep.*$", "", colnames(sub))
  col <- vl_palette_many_categ(length(unique(cditions)))[factor(cditions)]
  cor <- Hmisc::varclus(sub, similarity = "pearson", trans = "none")
  plot(cor, 
       hang = -1, 
       names.arg = NULL, 
       xlab = NA, 
       las = 1, 
       ylab= "Pearson r (25kb bins)")
  points(seq(ncol(sub)), 
         rep(0, ncol(sub)), 
         col = col[cor$hclust$order], 
         pch = 19)
}
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
bins <- vl_binBSgenome("dm6", 
                       bins_width = 25e3, 
                       steps_width = 25e3, 
                       restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY"))
dat <- meta[, .(score= vl_bw_coverage(bins, bw)), .(ChIP, cdition, rep)]
dat[, idx:= seq(.N), setdiff(names(dat), "score")]
mat <- dcast(dat,
             idx~ChIP+cdition+rep,
             value.var = "score")
mat <- as.matrix(mat, 1)

# Compute
cditions <- gsub("_rep.*$", "", colnames(mat))
col <- vl_palette_many_categ(length(unique(cditions)))[factor(cditions)]
cor <- Hmisc::varclus(mat, similarity = "pearson", trans = "none")

# Plot
pdf("pdf/cutnrun/PCC_cutnrun_replicates.pdf",
    width = 12,
    height = 5.5)
par(las= 1)
plot(cor, 
     hang = -1, 
     names.arg = NULL, 
     xlab = NA, 
     las = 1, 
     ylab= "Pearson r (25kb bins)")
points(seq(ncol(mat)), 
       rep(0, ncol(mat)), 
       col = col[cor$hclust$order], 
       pch = 19)
dev.off()
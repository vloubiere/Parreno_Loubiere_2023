setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(cdition, system, FC_file)]
dat[, cdition:= factor(cdition, 
                       levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dmat <- dcast(dat, FBgn+cdition~system, value.var = "diff")

pdf("pdf/Figures/Overlap_GFP_noGFP.pdf",5,5)
par(las= 2,
    mar= c(9,8,2.5,4.5))
dmat[, {
  .t <- table(GFP, noGFP)
  .t <- matrix(.t, ncol = ncol(.t), dimnames = dimnames(.t))
  rownames(.t) <- paste0(rownames(.t), "_GFP")
  colnames(.t) <- paste0(colnames(.t), "_noGFP")
  chi <- chisq.test(.t)
  conf <- chi$residuals
  vl_heatmap(matrix(conf, ncol = ncol(conf), dimnames = dimnames(conf)), 
             cluster_rows= F, 
             cluster_cols= F,
             legend_title = "OR (log2)", 
             display_numbers = T,
             display_numbers_matrix = .t,
             col= c("royalblue2", "yellow"),
             auto_margins= F,
             main= paste0("Chi. square pval= ", formatC(chi$p.value, format= "e")))
  print("done")
}, cdition]
dev.off()
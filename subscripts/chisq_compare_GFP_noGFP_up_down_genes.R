setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import data ----
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file), fread(FC_file), .(cdition, system, FC_file)]
dat[, cdition:= factor(cdition, 
                       levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dmat <- dcast(dat,
              FBgn+cdition~system,
              value.var = "diff")

# Plot ----
pdf("pdf/RNA_chisq_overlap_GFP_noGFP_systems.pdf", 12, 3)
vl_par(mfrow= c(1, 4))
dmat[, {
  .t <- table(GFP, noGFP)
  .t <- matrix(.t, ncol = ncol(.t), dimnames = dimnames(.t))
  rownames(.t) <- paste0(rownames(.t), "_GFP")
  colnames(.t) <- paste0(colnames(.t), "_noGFP")
  chi <- chisq.test(.t)
  conf <- chi$stdres
  vl_heatmap(matrix(conf, ncol = ncol(conf), dimnames = dimnames(conf)), 
             cluster.rows= F, 
             cluster.cols= F,
             legend.title = "Standardized\nresiduals\n", 
             display.numbers = T,
             display.numbers.matrix = .t,
             col= viridis::viridis(10),
             main= paste0(cdition, " (Chi.sq. pval= ", formatC(chi$p.value, format= "e"), ")"))
  print("done")
}, cdition]
dev.off()
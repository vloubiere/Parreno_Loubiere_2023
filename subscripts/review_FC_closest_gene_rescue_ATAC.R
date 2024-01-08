setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(glmnet)

# Import peaks ----
dat <- readRDS("Rdata/final_ATAC_rescue_table.rds")
peaks <- vl_resizeBed(dat, "center", 0, 0)
peaks <- peaks[diff_zfh1_ph=="unaffected"]
peaks[, class:= fcase(diff_zfh1_ph_vs_gfp_ph=="up", "Increased upon zfh1-KD",
                      diff_zfh1_ph_vs_gfp_ph=="unaffected", "Unaffected upon zfh1-KD",
                      diff_zfh1_ph_vs_gfp_ph=="down", "Decreased upon zfh1-KD")]
peaks[, class:= factor(class, 
                       c("Increased upon zfh1-KD",
                         "Unaffected upon zfh1-KD",
                         "Decreased upon zfh1-KD"))]

# Import genes ----
genes <- readRDS("Rdata/final_gene_features_table.rds")
TSS <- vl_resizeBed(genes, "start", 0, 0)
cl <- vl_closestBed(peaks, TSS)
cl <- cl[abs(dist)<25e3]

# Plot
Cc <- c("tomato", "lightgrey", "cornflowerblue")
Cc <- adjustcolor(Cc, .6)

pdf("pdf/review_boxplot_FC_rescued_ATAC_peaks_closest_genes.pdf", 4.4, 3)
vl_par(mai= c(.9, 2.9, .9, .9),
       mgp= c(1, .25, 0),
       lwd= .5)
cl[, {
  vl_boxplot(log2FoldChange_GFP_PH.b~class,
             xlab= "Closest gene fold change (log2)\n gfp+ph-KD vs gfp+w-KD",
             compute.pval= list(c(1,2), c(2,3)),
             col= Cc,
             tilt.names= T,
             horizontal= T,
             xaxt= "n")
  axis(1, padj= -1.25)
  abline(v= 0,
         lty= "11")
}]
dev.off()
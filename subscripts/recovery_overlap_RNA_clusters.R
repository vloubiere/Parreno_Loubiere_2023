setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/final_gene_features_table.txt")
dat[, cl:= as.character(cl)]
dat[!is.na(cl), cl:= paste0("cluster ", cl)]
dat[is.na(cl), cl:= "None"]

###########################################
# PLOT
###########################################
pdf("pdf/recovery_overlap_RNA_clusters.pdf", 
    width = 5,
    height = 5)
par(las= 1,
    mar= c(8,7,8,4),
    tcl= -0.2,
    mgp= c(2,0.5,0))
dat[!is.na(recovery), {
  .t <- table(recovery, cl, useNA = "ifany")
  mat <- matrix(.t, ncol= ncol(.t), dimnames = dimnames(.t))
  vl_heatmap(mat, 
             cluster_rows= F,
             cluster_cols= F, 
             display_numbers= T, 
             show_legend= F)
}]
dev.off()
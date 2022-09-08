setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/final_gene_features_table.txt")
dat[, cl:= as.character(cl)]
dat[!is.na(cl), cl:= paste0("cluster ", cl)]
dat[is.na(cl), cl:= "None"]
dat <- dat[!is.na(recovery)]

###########################################
# PLOT
###########################################
# Ovelap with RNA clusters
pdf("pdf/Figures/characterization_recovery_groups.pdf", 
    width = 5,
    height = 5)
par(las= 1,
    mar= c(8,7,8,4))
dat[!is.na(recovery), {
  .t <- table(recovery, cl, useNA = "ifany")
  mat <- matrix(.t, ncol= ncol(.t), dimnames = dimnames(.t))
  vl_heatmap(mat, 
             cluster_rows= F,
             cluster_cols= F, 
             display_numbers= T, 
             show_legend= F)
}]

# log2FC per cdition
.m <- melt(dat, id.vars = "recovery", measure.vars = patterns("log2FoldChange"))
.m[, variable:= gsub("^log2FoldChange_", "", variable)]
par(mar= c(7,7,3,5))
test <- vl_boxplot(value~recovery+variable, 
           .m, ylab= "log2FoldChange",
           col= c("rosybrown1", "palegreen3"),
           compute_pval= list(c(1,2), c(3,4), c(5,6), c(7,8)),
           xaxt= "n",
           tilt.names= T,
           at= rep(seq(1,7,2), each= 2)+c(0.2,0.8))
axis(1, 
     at= c(1.5, 3.5, 5.5, 7.5),
     c("PH18", "PH29", "PHD11", "PHD9"))
legend("topright",
       fill= c("palegreen3", "rosybrown1"),
       legend= c("Recovery", "No recovery"),
       bty= "n")
dev.off()
setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/final_gene_features_table.txt")

###########################################
# PLOT
###########################################
pdf("pdf/recovery_overlap_K27_K118_domains.pdf", 
    width = 5,
    height = 5)
par(las= 1,
    mar= c(8,14,8,6),
    tcl= -0.2,
    mgp= c(2,0.5,0))
dat[!is.na(recovery), {
  vl_upset_plot(list(K27me3_bound= .SD[(K27me3_bound), FBgn],
                     K118Ub_bound= .SD[(K118Ub_bound), FBgn]))
  title(main= recovery)
}, recovery]
dev.off()
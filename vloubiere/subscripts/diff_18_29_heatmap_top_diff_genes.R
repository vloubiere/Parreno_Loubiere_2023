setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

#------------------#
# Heatmap selected genes
#------------------#
FC <- readRDS("Rdata/final_FC_table.rds")
FC <- FC[grepl("dose", cdition)]

dat <- readRDS("Rdata/clustering_dose_18_29_transcriptomes.rds")
dat <- unique(dat[, .(cl= rcl, diff= range(value)[2]-range(value)[1]), .(FBgn= row)])
dat <- dat[order(diff, decreasing = T)][1:100]
mat <- dcast(FC[FBgn %in% dat$FBgn], 
             symbol~cdition, 
             value.var = "log2FoldChange")
mat <- as.matrix(mat, 1)

pdf("pdf/comparison_ph18_ph29_DOSE/heatmap_top_diff_genes.pdf", height = 20)
vl_heatmap(mat, 
           cluster_cols = F, 
           breaks = c(-8, 0, 8), 
           display_numbers = T, 
           legend_title = "log2FoldChange")
dev.off()
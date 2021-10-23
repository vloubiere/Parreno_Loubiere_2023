setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

FC <- readRDS("Rdata/final_FC_table.rds")
FC <- FC[grepl("dose", cdition)]
FC <- dcast(FC, 
            symbol~cdition, 
            value.var = "log2FoldChange")
setkeyv(FC, "symbol")

#-------------------------#
# Motif enrichment kmeans clusters
#-------------------------#
if(!file.exists("db/motif_enrichment_tables/dose_diff_18_29.rds"))
{
  cl_diff <- readRDS("Rdata/clustering_dose_18_29_transcriptomes.rds")
  cl_diff <- unique(cl_diff[, .(FBgn= row, rcl)])
  if(!exists("mot"))
    mot <- readRDS("Rdata/ED_REs_motif_counts_final.rds")
  mot[cl_diff, cl:= i.rcl, on= "FBgn"]
  mot[is.na(cl), cl:= 0]
  enr_all <- vl_motif_cl_enrich(mot, 
                                cl_column = "cl", 
                                bg = 0)
  symbols <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
  enr_all[symbols, FBgn:= i.gene_id, on= "motif_name==gene_symbol"]
  saveRDS(enr_all, "db/motif_enrichment_tables/dose_diff_18_29.rds")
}
enr_all <- readRDS("db/motif_enrichment_tables/dose_diff_18_29.rds")
enr_all <- enr_all[motif %in% enr_all[order(padj), motif[1], motif_name]$V1]



pdf("pdf/comparison_ph18_ph29_DOSE/motif_enrichment_clusters.pdf", 
    width = 5, 
    height = 10)
layout(matrix(1:2, ncol= 2), 
       widths = c(1,0.6))
pl <- vl_motif_cl_enrich_plot_only(enr_all, 
                                   padj_cutoff = 0.01, 
                                   N_top = 8)
mat <- FC[unique(pl[order(y, decreasing = T), .(y, motif_name)]$motif_name)]
mat <- as.matrix(mat, 1)
par(mar= c(5, 4.25, 4.5, 6),
    cex= 0.5)
vl_heatmap(mat, 
           cluster_rows = F, 
           auto_margins = F,  
           legend_title = "log2FC", 
           display_numbers = T, 
           breaks = c(-5, 0, 5))
box(lwd= 0.25)
dev.off()

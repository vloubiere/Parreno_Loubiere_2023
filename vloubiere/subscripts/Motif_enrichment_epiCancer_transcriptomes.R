setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

dat <- readRDS("Rdata/clustering_epiCancer_transcriptomes.rds")
FC <- readRDS("Rdata/final_FC_table.rds")
FC <- FC[grepl("TS", cdition)]
FC <- FC[!grepl("18", cdition)]
FC <- dcast(FC, 
            symbol~cdition, 
            value.var = "log2FoldChange")
setkeyv(FC, "symbol")

#------------------------------------#
# Compute motif enrichment
#------------------------------------#
# All REs
if(!file.exists("db/motif_enrichment_tables/epiCancer_clusters_all.rds"))
{
  mot <- readRDS("Rdata/ED_REs_motif_counts_final.rds")
  mot <- mot[dat, cl:= i.rcl, on= "FBgn"]
  mot[is.na(cl), cl:=0]
  res <- vl_motif_cl_enrich(obj = mot, 
                            cl_column = "cl", 
                            bg = 0)
  saveRDS(res,
          "db/motif_enrichment_tables/epiCancer_clusters_all.rds")
}

enr_all <- readRDS("db/motif_enrichment_tables/epiCancer_clusters_all.rds")
enr_all <- enr_all[motif %in% enr_all[order(padj), motif[1], motif_name]$V1]

pdf("pdf/clustering/motif_enrichment_epiCancer_clusters_all.pdf", 
    height = 12, 
    width = 7)
layout(matrix(1:2, ncol=2), 
       widths = c(1, 0.5))
pl <- vl_motif_cl_enrich_plot_only(enr_all,
                                   padj_cutoff = 0.01, 
                                   N_top = 10)
mat <- FC[unlist(tstrsplit(unique(pl[order(y, decreasing = T), .(y, motif_name)]$motif_name), "/", keep=1))]
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

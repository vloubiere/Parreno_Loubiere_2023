setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

#------------------------------------------------#
# Import data
#------------------------------------------------#
dat <- fread("Rdata/cl2_cl5_RECOVERY_genes.txt")
ATAC <- readRDS("Rdata/motif_counts_ATAC_sites.rds")
dat <- ATAC[dat, on= "FBgn"]

#------------------------------------------------#
# Plot
#------------------------------------------------#
pdf("pdf/Figures/motif_enrich_cl2_vs_cl5.pdf", 
    width = 12, 
    height = 8)
par(mar= c(5,10,2,8),
    mfrow= c(1,3))
dat[, {
  # Motif enrichment
  enr <- vl_motif_enrich(counts= as.matrix(.SD[(RECOVERY)]),
                         control_counts = as.matrix(.SD[!(RECOVERY)]),
                         plot= F)
  enr[, variable:= gsub("_mot$", "", variable)]
  # Collapse per motif cluster
  enr[vl_Dmel_motifs_DB_full, motif_cluster:= Motif_cluster_name, on= "variable==uniqName_noSpecialChar"]
  enr[, sel:= seq(.N)==which.min(padj), .(log2OR>0, motif_cluster)]
  enr <- enr[(sel)]
  # plot
  plot(enr, padj_cutoff = 0.05)
  abline(v= 0, lty= "11")
  title(main= paste0("Recov. (cl5) vs not (cl2) ->", group))
  print(".")
}, group, .SDcols= patterns("_mot$")]
dev.off()

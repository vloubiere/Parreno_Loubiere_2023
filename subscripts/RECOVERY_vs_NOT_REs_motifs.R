setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

#------------------------------------------------#
# Import data
#------------------------------------------------#
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
dat <- dat[, .(FBgn, PHD9_RECOVERY, PHD11_RECOVERY)]
ATAC <- readRDS("Rdata/motif_counts_ATAC_sites.rds")
ATAC <- ATAC[, !c("coor", "symbol")]
dat <- ATAC[dat, on= "FBgn"]
dat <- melt(dat, 
            measure.vars = c("PHD9_RECOVERY", "PHD11_RECOVERY"), 
            value.name = "class")
dat <- dat[!is.na(class)]

#------------------------------------------------#
# Plot
#------------------------------------------------#
pdf("pdf/Figures/motif_enrich_recovery_vs_no.pdf", 
    width = 9, 
    height = 8)
par(mar= c(5,10,2,8),
    mfrow= c(1,2))
dat[, {
  # Motif enrichment
  enr <- vl_motif_enrich(counts= as.matrix(.SD[(class)]),
                         control_counts = as.matrix(.SD[!(class)]),
                         plot= F)
  enr[, variable:= gsub("_mot$", "", variable)]
  # Collapse per motif cluster
  enr[vl_Dmel_motifs_DB_full, motif_cluster:= Motif_cluster_name, on= "variable==uniqName_noSpecialChar"]
  enr[, sel:= seq(.N)==which.min(padj), .(log2OR>0, motif_cluster)]
  enr <- enr[(sel)]
  # plot
  plot(enr, padj_cutoff = 0.01)
  abline(v= 0, lty= "11")
  title(main= paste(gsub("_RECOVERY$", "", variable), "recovery vs no recovery", group))
  print(".")
}, .(group, variable), .SDcols= patterns("_mot$")]
dev.off()
setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

#------------------------------------------------#
# Import data
#------------------------------------------------#
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
dat <- dat[, .(FBgn, PHD9_RECOVERY, PHD11_RECOVERY)]
RE <- fread("Rdata/final_RE_motifs_table.txt")
dat <- RE[dat, on= "FBgn"]
dat <- melt(dat,
            measure.vars = c("PHD9_RECOVERY", "PHD11_RECOVERY"), 
            value.name = "class")
dat <- dat[!is.na(class)]

#------------------------------------------------#
# Plot
#------------------------------------------------#
pdf("pdf/Figures/motif_enrich_recovery_vs_no.pdf", 
    width = 9, 
    height = 4)
par(mar= c(5,10,2,8),
    mfrow= c(1,3))
dat[, {
  # Motif enrichment
  enr <- vl_motif_enrich(counts= as.matrix(.SD[(class)]),
                         control_counts = as.matrix(.SD[!(class)]),
                         collapse_clusters = vl_Dmel_motifs_DB_full[, .(motif= paste0(motif, "_mot"), motif_cluster)],
                         plot= F)
  # plot
  plot(enr, padj_cutoff = 0.01)
  abline(v= 0, lty= "11")
  title(main= paste(gsub("_RECOVERY$", "", variable), "recovery vs no recovery", group))
  print(".")
}, keyby= .(variable, group), .SDcols= patterns("_mot$")]
dev.off()
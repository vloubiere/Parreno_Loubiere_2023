setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]
dat <- vl_resizeBed(dat, "start", 750, 250)
sel <- vl_Dmel_motifs_DB_full[Dmel %in% c("pho", "Spps", "Stat92E") | motif_cluster=="STAT6", motif]
counts <- vl_motif_counts(dat[, .(seqnames, start, end)], 
                          sel= sel,
                          genome = "dm6",
                          p.cutoff = 1e-03)
dat <- cbind(dat, counts)
pl <- melt(dat, id.vars = "recovery", measure.vars = sel)
pl[, recovery:= factor(recovery, c("Recovery", "noRecovery"))]
setorderv(pl, "recovery")
pl[, col:= c("palegreen3", "rosybrown1")[.GRP], recovery]

pdf("pdf/recovery_promoters_motif_counts.pdf", 4, 4)
par(las= 1,
    tcl= -0.2,
    mgp= c(2,0.5,0))
pl[, {
  plot(NA,
       xlim= c(0, max(value)),
       ylim= c(0, max(table(value))),
       main= variable,
       xlab= "Motif counts",
       ylab= "Promoter counts")
  br <- seq(0, max(value))
  .SD[, {
    hist(value,
         col= adjustcolor(col[1], 0.5),
         add= T, 
         breaks= br)
    print("")
  }, .(recovery, col)]
  legend("topright",
         legend= c("Recovery (68)",
                   "noRecovery (61)"),
         fill= adjustcolor(c("palegreen3", "rosybrown1"), 0.5),
         bty= "n")
  print("")
}, variable]
dev.off()

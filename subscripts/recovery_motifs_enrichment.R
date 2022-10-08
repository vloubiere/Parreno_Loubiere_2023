setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]
motifs <- fread("Rdata/final_RE_motifs_table.txt")
motifs <- motifs[between(dist, -5000, 0)] # REs distance cutoff
motifs <- motifs[dat, on="FBgn"]
mot <- names(motifs)[names(motifs) %in% vl_Dmel_motifs_DB_full$motif]
motifs <- motifs[, lapply(.SD, sum), .(group, recovery, FBgn), .SDcols= mot]

# Compute motif enrichments
groups <- split(motifs, motifs$group)
enr <- lapply(groups, function(x) {
  .c <- vl_motif_enrich(counts = x[recovery=="Recovery", ..mot],
                  control_counts = x[recovery=="noRecovery", ..mot],
                  collapse_clusters = vl_Dmel_motifs_DB_full[mot, motif_cluster, on= "motif"],
                  plot= F)
  .c[vl_Dmel_motifs_DB_full, TF:= i.Dmel, on= "variable==motif_cluster", mult= "first"]
  .c[!is.na(TF), variable:= ifelse(grepl(TF, variable), variable, paste0(variable, " -> ", TF)), .(variable, TF)]
  return(.c)
})

#############################
# PLOT
#############################
pdf("pdf/recovery_motifs_enrichment.pdf", 6, 4)
par(mar= c(4,40,3,7),
    mgp= c(2,0.5,0),
    las= 1,
    cex= 0.5,
    tcl= -0.2)
for(i in seq(enr))
{
  pl <- plot(enr[[i]], 
             padj_cutoff= 0.05, 
             top_enrich= 20, order= "log2OR")
  vl_add_motifs(pl)
  title(main= paste0("Motif enrichment ", names(enr)[i], " Recovery/No recovery"))
}
dev.off()

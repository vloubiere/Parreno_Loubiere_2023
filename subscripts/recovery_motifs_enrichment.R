setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(cl)]
dat <- dat[!is.na(recovery)]
motifs <- fread("Rdata/final_RE_motifs_table.txt")
motifs <- dat[motifs, on="FBgn", nomatch= NULL]

# Compute motif enrichments
groups <- split(motifs, motifs$group)
mot <- names(motifs)[names(motifs) %in% vl_Dmel_motifs_DB_full$motif]
enr <- lapply(groups, function(x) {
  vl_motif_enrich(counts = as.matrix(x[recovery=="Recovery", ..mot]),
                  control_counts = as.matrix(x[recovery!="Recovery", ..mot]),
                  vl_Dmel_motifs_DB_full[mot, Dmel, on= "motif"],
                  plot= F)
})

#############################
# PLOT
#############################
pdf("pdf/recovery_motifs_enrichment.pdf", 3, 5)
par(mar= c(4,10,3,7),
    mgp= c(2,0.5,0),
    las= 1,
    cex= 0.5,
    tcl= -0.2)
for(i in seq(enr))
{
  plot(enr[[i]],
       top_enrich= 40)
  title(main= paste0("Motif enrichment ", names(enr)[i], " Recovery/No recovery"))
}
dev.off()

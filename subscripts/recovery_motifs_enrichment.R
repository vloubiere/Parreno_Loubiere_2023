setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]

# Compute motif enrichments
sel <- names(motifs)[names(motifs) %in% vl_Dmel_motifs_DB_full$motif]
# motifs <- motifs[group=="TSS"]
enr <- motifs[, {
  # vl_motif_cl_enrich(counts_matrix = .SD, 
  #                    cl_IDs = cl,
  #                    control_cl = "0",
  #                    collapse_clusters = vl_Dmel_motifs_DB_full[sel, motif_cluster, on= "motif"],
  #                    plot= F)
  vl_motif_enrich(.SD[cl=="Recovery"],
                  .SD[cl!="Recovery"], 
                  collapse_clusters = vl_Dmel_motifs_DB_full[sel, motif_cluster, on= "motif"],
                  plot= F)
}, group, .SDcols= sel]

#############################
# PLOT
#############################
pdf("pdf/recovery_motifs_enrichment.pdf", 6, 15)
par(mar= c(4,40,3,7),
    mgp= c(2,0.5,0),
    las= 1,
    cex= 0.5,
    tcl= -0.2)
enr[padj<0.05 & log2OR>0, {
  .c <- .SD
  setattr(.c, "class", c("vl_enr", "data.table", "data.frame"))
  plot(.c, 
       padj_cutoff= 0.05)
  ""
}, group]
dev.off()

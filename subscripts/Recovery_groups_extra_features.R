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
dat[, col:= ifelse(recovery=="Recovery", "palegreen3", "rosybrown1")]
motifs <- fread("Rdata/final_RE_motifs_table.txt")
motifs <- dat[motifs, on="FBgn", nomatch= NULL]

# Compute gene network
cols <- c("log2FoldChange_PH29", "log2FoldChange_PHD9", "log2FoldChange_PHD11")
sizes <- apply(dat[, ..cols], 1, function(x) max(abs(x), na.rm=T))
sizes[sizes>quantile(sizes, 0.90)] <- quantile(sizes, 0.90)
sizes <- log2(sizes+2)
net <- vl_STRING_interaction(symbols = dat$symbol,
                             species = "Dm",
                             col= dat$col,
                             size = sizes*2.5,
                             cex.label = sizes/4,
                             plot= F)

# Compute GO enrich
GO <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat$recovery), 
                   species = "Dm", 
                   plot= F)

# Compute motif enrichments
groups <- split(motifs, motifs$group)
mot <- names(motifs)[names(motifs) %in% vl_Dmel_motifs_DB_full$motif]
enr <- lapply(groups, function(x) {
  vl_motif_enrich(counts = as.matrix(x[recovery=="Recovery", ..mot]),
                  control_counts = as.matrix(x[recovery!="Recovery", ..mot]),
                  plot= F)
})

saveRDS(list(net= net, 
             GO= GO,
             enr= enr), 
        "Rdata/recovery_group_features.rds")

#############################
# PLOT
#############################
pdf("pdf/Figures/recovery_groups_features.pdf", 7, 5)
par(mar= c(3,31,3,7),
    las= 1,
    cex= 0.8)
# GOs
plot(GO, 
     padj_cutoff = 1e-5, 
     top_enrich = 20, 
     cex.balloons= 0.5)

# Network
par(mar= c(1,1,1,1))
set.seed(1)
plot(net)
legend("topleft",
       fill= c("palegreen3", "rosybrown1"),
       legend= c("Recovery", "No recovery"),
       bty= "n")
# Motifs
par(mar= c(4,25,3,22),
    cex= 0.5)
for(i in seq(enr))
{
  plot(enr[[i]],
       top_enrich= 40)
  title(main= paste0("Motif enrichment ", names(enr)[i], " Recovery/No recovery"))
}
dev.off()
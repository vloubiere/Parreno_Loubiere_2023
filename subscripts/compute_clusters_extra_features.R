setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(cl)]
motifs <- fread("Rdata/final_RE_motifs_table.txt")
motifs <- motifs[between(dist, -5000, 0)]
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
GO_all <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat$cl),
                       geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn,
                       species = "Dm", 
                       plot= F)
GO_PRC1 <- vl_GO_enrich(geneIDs = split(dat$FBgn, 
                              dat[, .(cl, ifelse(PRC1_bound, "PRC1+", "PRC1-"))]),
                        geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn,
                        species = "Dm", 
                        plot= F)

# Compute motif enrichments
groups <- split(motifs, motifs$group)
mot <- names(motifs)[names(motifs) %in% vl_Dmel_motifs_DB_full$motif]
enr <- lapply(groups, function(x) {
  vl_motif_cl_enrich(counts_matrix = as.matrix(x[, ..mot]),
                     collapse_clusters = vl_Dmel_motifs_DB_full[mot, Dmel, on= "motif"],
                     cl_IDs = x[, paste0(cl, ifelse(PRC1_bound, " PRC1+", " PRC1-"))],
                     plot= F)
})
saveRDS(list(net= net, 
             GO_all= GO_all,
             GO_PRC1= GO_PRC1,
             enr= enr), 
        "Rdata/clustering_RNA_features.rds")

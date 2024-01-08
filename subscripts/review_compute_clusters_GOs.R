setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data
dat <- readRDS("Rdata/final_gene_features_table.rds")
universe <- dat$FBgn
dat <- dat[cl!="Unaffected"]

# GO enrichment ----
GO_PcG <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat[, paste(cl, ifelse(PcG_bound, "PcG+", "PcG-"))]),
                       geneUniverse.IDs = universe,
                       species = "Dm",
                       select = c("BP", "MF"),
                       plot= F)

# Save
saveRDS(GO_PcG,
        "Rdata/clustering_RNA_GOs.rds")
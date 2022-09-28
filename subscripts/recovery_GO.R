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

# Compute GO enrich
GO <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat$recovery), 
                   geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn,
                   species = "Dm",
                   plot= F)

#############################
# PLOT
#############################
pdf("pdf/recovery_GOs.pdf", 4, 5)
par(mar= c(3,31,3,7),
    las= 1,
    cex= 0.5)
plot(GO, 
     padj_cutoff = 0.05, 
     top_enrich = 40, 
     cex.balloons= 0.4)
dev.off()

# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[is.na(cl), cl:= 0]

#------------------------------#
# Compute motif enrichments
#------------------------------#
TSS <- vl_resizeBed(dat[, .(seqnames, start, end, strand)], 
                    center = "start", 
                    upstream = 200, 
                    downstream = 50)
counts <- vl_motif_counts(TSS, 
                          genome = "dm6", 
                          sel= vl_Dmel_motifs_DB_full[collection %in% c("bergman", "flyfactorsurvey", "jaspar") & !is.na(FBgn), motif_ID])
counts_list <- split(counts, dat[, paste(cl, ifelse(PRC1_bound|K27me3_bound|K118Ub_bound, "PcG+", "PcG-"))])
enr <- vl_motif_cl_enrich(counts_list = counts_list, control_cl = c("0 PcG-", "0 PcG+"))

# Collapse Dmel ID
enr[vl_Dmel_motifs_DB_full, name:= i.Dmel, on= "variable==motif_ID"]
enr <- enr[variable %in% enr[, .SD[log2OR>0][which.min(padj), variable], name]$V1]

#------------------------------#
# Compute network and GO
#------------------------------#
# Restrict to clustered genes
dat <- dat[cl!=0]
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
                             plot= F, 
                             version = "11.0")

# Compute GO enrich
GO_all <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat$cl),
                       geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn,
                       species = "Dm", 
                       plot= F)
GO_PcG <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat[, paste(cl, ifelse(PRC1_bound|K27me3_bound|K118Ub_bound, "PcG+", "PcG-"))]),
                       geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn,
                       species = "Dm", 
                       plot= F)

# Save
saveRDS(list(net= net, 
             GO_all= GO_all,
             GO_PcG= GO_PcG,
             enr= enr), 
        "Rdata/clustering_RNA_features.rds")

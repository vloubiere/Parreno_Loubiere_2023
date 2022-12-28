# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[is.na(cl), recovery:= 0]
dat <- dat[!is.na(recovery)]

#------------------------------#
# Compute motif enrichments
#------------------------------#
TSS <- vl_resizeBed(dat[, .(seqnames, start, end, strand)], 
                    center = "start", 
                    upstream = 200, 
                    downstream = 50)
sel <- vl_Dmel_motifs_DB_full[!is.na(FBgn), motif_ID]
counts <- vl_motif_counts(TSS, 
                          genome = "dm6", 
                          sel= sel)
counts_list <- split(counts, dat$recovery)
enr <- vl_motif_cl_enrich(counts_list = counts_list, control_cl = "0")

# Collapse Dmel ID
enr[vl_Dmel_motifs_DB_full, name:= i.Dmel, on= "variable==motif_ID"]
enr <- enr[variable %in% enr[, .SD[log2OR>0][which.min(padj), variable], name]$V1]

pdf("pdf/Figure_3_motifs_enrichment.pdf", 6.5, 8)
par(mar= c(4,20,3,7),
    mgp= c(2,0.5,0),
    las= 1,
    tcl= -0.2)
plot(enr, 
     padj_cutoff= 0.05, 
     top_enrich= 30)
dev.off()


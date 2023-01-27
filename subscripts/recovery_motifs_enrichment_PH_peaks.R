# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt", 
             sel= c("seqnames", "start", "end", "strand", "recovery", "FBgn"))
dat[is.na(recovery), recovery:= "control"]

#------------------------------#
# Assign PH peaks to genes
#------------------------------#
PH <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
PH[, start:= start+peak]
PH[, end:= start]
TSS <- vl_resizeBed(dat, 
                    center = "start", 
                    upstream = 0, 
                    downstream = 0)
PH <- vl_closestBed(PH, TSS)
PH <- PH[, .(seqnames, start, end, strand, name, FBgn= FBgn.b, recovery= recovery.b)]
PH <- PH[!is.na(recovery)]

# Remove peaks that are too distal or that could not be unambiguously assigned
genes <- vl_resizeBed(dat, 
                      center = "start", 
                      upstream = 5000, 
                      downstream = dat[, end-start])
PH[genes, c("gene_start", "gene_end"):= .(i.start, i.end), on= "FBgn"]
PH <- PH[between(start, gene_start, gene_end)]
PH <- PH[, check:= length(unique(recovery))==1, name][(check), !"check"]

# Resize
PH <- vl_resizeBed(PH, "center", 125, 125)

#-----------------------------------#
# Motif enrichment
#-----------------------------------#
sel <- vl_Dmel_motifs_DB_full[!is.na(FBgn), motif_ID]
counts <- vl_motif_counts(PH, 
                          genome = "dm6", 
                          sel= sel)
counts_list <- split(counts, PH$recovery)
enr <- vl_motif_cl_enrich(counts_list = counts_list, control_cl = "control")
enr[vl_Dmel_motifs_DB_full, name:= i.Dmel, on= "variable==motif_ID"]

pdf("pdf/Figure_3_motifs_enrichment_PH_peaks.pdf", 6.5, 9)
par(mar= c(4,24,3,7),
    mgp= c(2,0.5,0),
    las= 1,
    tcl= -0.2)
plot(enr[variable %in% enr[, .SD[log2OR>0][which.min(padj), variable], name]$V1], 
     padj_cutoff= 0.01, 
     top_enrich= 30,
     cex.balloons= 0.8)
dev.off()


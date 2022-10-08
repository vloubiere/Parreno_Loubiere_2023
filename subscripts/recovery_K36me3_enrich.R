setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
gf <- fread("Rdata/final_gene_features_table.txt")

# Format
.m <- melt(gf, 
            id.vars = c("FBgn", "recovery"), 
            measure.vars = patterns("FPKM"= "^FPKM_PH"), 
           value.name = "FPKM")
.m[, cdition:= tstrsplit(variable, "_", keep= 2)]
.m$variable <- NULL
setorderv(.m, c("cdition", "FPKM"))
.m[, idx:= .I, cdition]
.m[, ctl:= idx+1, cdition]
while(any("noRecovery" %in% .m[.m[recovery=="noRecovery", ctl], recovery]))
  .m[recovery=="noRecovery", ctl:= ctl+1, cdition]

# Select activity matched controls
dat <- rbind(.m[recovery=="noRecovery", FBgn:cdition],
             .m[.m[recovery=="noRecovery", ctl], FBgn:cdition])
dat[recovery!="noRecovery" | is.na(recovery), recovery:= "control"]

# Add K36me3
K36 <- melt(gf, 
            id.vars = c("FBgn", "recovery"), 
            measure.vars = patterns("^H3K36me3"), 
            value.name = "K36me3")
K36[, cdition:= tstrsplit(variable, "_", keep= 2)]
dat[K36, K36me3:= K36me3, on= c("FBgn", "cdition")]

# Add Tracks
dat[, track:= list.files("db/bw/cutnrun/", 
                         paste0("H3K36me3_", cdition, "_merge"), 
                         full.names = T), cdition]

# Add TTS
gtf <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)[type=="gene"]
TTS <- vl_resizeBed(gtf, "end", 0, 0)
TTS <- TTS[, .(gene_id, seqnames, start, end, strand)]
dat <- merge(dat, TTS, by.x= "FBgn", by.y= "gene_id")
dat[, recovery:= factor(recovery, c("control", "noRecovery"))]

# PLOT
Cc <- c("lightgrey", "rosybrown1")

pdf("pdf/recovery_H3K36me3_noRecovery_genes.pdf", 9, 5)
layout(matrix(1:12, 4, 3, byrow = T), 
       widths = c(0.35,0.35,1))
par(mar= c(5,3,2,1),
    mgp= c(1.5,0.5,0),
    tcl= -0.2,
    las= 1)
dat[, {
  vl_boxplot(FPKM~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2)),
             main= cdition,
             col= Cc,
             ylab= "FPKM")
  vl_boxplot(K36me3~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2)),
             main= cdition,
             col= Cc,
             ylab= "Enrichment")
  vl_bw_average_track(bed = .SD, 
                      tracks = track, 
                      upstream = 10000,
                      downstream = 10000,
                      set_IDs = recovery, 
                      stranded = T, 
                      center_label = "TTS", 
                      col= Cc)
  print("")
}, .(cdition, track), .SDcols= c("seqnames", "start", "end", "strand")]
dev.off()
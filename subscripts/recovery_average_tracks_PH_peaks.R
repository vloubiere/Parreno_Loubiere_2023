setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(rtracklayer)
require(GenomicRanges)
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]
dat <- melt(dat, 
            id.vars = c("FBgn", "recovery"), 
            measure.vars = patterns("_body$|_prom$|_TTS$"))
dat[, variable:= gsub("_body$|_prom$|_TTS$", "", variable)]
dat[, track:= list.files("db/bw/cutnrun/", 
                         paste0(variable, "_merge"), 
                         full.names = T), variable]
dat[is.na(track), track := list.files("db/bw/SA_2020/", 
                                      paste0(variable, "_merge"), 
                                      full.names = T), variable]
dat[variable=="ATAC", track:= "db/bw/ATAC/ATAC_merged.bw"]
# TSS coor
peaks <- fread("Rdata/final_RE_motifs_table.txt")[group=="PH", .(coor, FBgn)]
peaks[, c("seqnames", "start", "end"):= tstrsplit(coor, ":|-")]
dat <- merge(dat, peaks, by= "FBgn")
# Factors
dat[, recovery:= factor(recovery, c("Recovery", "noRecovery"))]
dat[, c("variable", "cdition"):= tstrsplit(variable, "_")]
dat[variable=="ATAC", cdition:= "ED"]
dat[, cdition:= factor(cdition, c("PH18", "PH29", "PHD9", "PHD11", "ED", "CNSID"))]

#-----------------------------------------#
# Plot
#-----------------------------------------#
Cc <- c("palegreen3", "rosybrown1")

pdf("pdf/recovery_ph18_average_tracks_PH_peaks.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
dat[!(cdition %in% c("PH29", "PHD9", "PHD11")), {
  # Plot average track for each class
  .c <- unique(.SD)
  .q <- vl_bw_average_track(bed= .c, 
                            tracks= track,
                            set_IDs = .c$recovery,
                            upstream = 2500,
                            downstream = 2500,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            legend= F,
                            col = Cc)
  
  # Legend
  leg <- .c[, .(legend= paste0(recovery, " (", .N, " sites)")), keyby= recovery]
  legend("topright",
         legend = leg$legend,
         fill= Cc,
         bty= "n")
  title(main= paste(variable, cdition))
  
  # Boxplot
  vl_boxplot(value~recovery,
             col= Cc,
             ylab= "Enrichment", 
             compute_pval= list(c(1,2)),
             tilt.names= T)
  print(".")
}, .(variable, cdition, track), .SDcols= c("seqnames", "start", "end", "recovery")]
dev.off()

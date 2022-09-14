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
# TSS coor
gtf <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
TSS <- vl_resizeBed(gtf[type=="gene"], "start", 0, 0)
TSS <- TSS[, .(seqnames, start, end, strand, gene_id)]
dat <- merge(dat, TSS, by.x= "FBgn", by.y= "gene_id")
# Factors
dat[, recovery:= factor(recovery, c("Recovery", "noRecovery"))]
dat[, c("variable", "cdition"):= tstrsplit(variable, "_")]
dat[, cdition:= factor(cdition, c("PH18", "PH29", "PHD9", "PHD11", "ED"))]

#-----------------------------------------#
# Plot 1
#-----------------------------------------#
Cc <- c("palegreen3", "rosybrown1")

pdf("pdf/Figures/recovery_ph18_average_tracks.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
dat[cdition %in% c("PH18", "ED"), {
  # Plot average track for each class
  .c <- unique(.SD)
  .q <- vl_bw_average_track(bed= .c, 
                            tracks= track,
                            set_IDs = .c$recovery,
                            upstream = 2500,
                            downstream = 10000,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            legend= F,
                            col = Cc)
  
  # Legend
  leg <- .c[, .(legend= paste0(recovery, " (", .N, " genes)")), keyby= recovery]
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
}, .(variable, cdition, track), .SDcols= c("seqnames", "start", "end", "strand", "recovery")]
dev.off()

#-----------------------------------------#
# Plot 2
#-----------------------------------------#
Cc <- c("chartreuse3", "chocolate1", "darkorchid2", "brown2")

pdf("pdf/Figures/recovery_cditions_average_tracks.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
dat[cdition!="ED", {
  # Plot average track for each class
  .c <- unique(.SD)
  .c[, cdition:= factor(cdition, c("PH18", "PH29", "PHD9", "PHD11"))]
  setorderv(.c, "cdition")
  vl_bw_average_track(bed= .c, 
                      set_IDs = .c$cdition,
                      tracks= unique(track),
                      upstream = 2500,
                      downstream = 10000,
                      stranded = T,
                      center_label = "TSS", 
                      legend.cex = 0.6, 
                      col = Cc,
                      legend= F,
                      col.adj= c(0.3,0.7))
  # Legend
  legend("topright",
         legend = levels(.c$cdition),
         fill= Cc,
         bty= "n",
         cex= 0.7)
  title(main= paste(recovery, variable))
  
  # Boxplot
  vl_boxplot(value~cdition,
             .c,
             col= Cc,
             ylab= "Enrichment", 
             compute_pval= list(c(1,2), c(1,3), c(1,4)),
             tilt.names= T)
  print(".")
}, .(recovery, variable), .SDcols= c("seqnames", "start", "end", "strand", "cdition", "track")]
dev.off()
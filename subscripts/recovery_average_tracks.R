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
gtf <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
TSS <- vl_resizeBed(gtf[type=="gene"], "start", 0, 0)
TSS <- TSS[, .(seqnames, start, end, strand, gene_id)]
dat <- merge(dat, TSS, by.x= "FBgn", by.y= "gene_id")
# Factors
dat[, recovery:= factor(recovery, c("Recovery", "noRecovery"))]
dat[, c("variable", "cdition"):= tstrsplit(variable, "_")]
dat[variable=="ATAC", cdition:= "ED"]

#-----------------------------------------#
# Plot 1
#-----------------------------------------#
Cc <- c("palegreen3", "rosybrown1")

pdf("pdf/recovery_ph18_average_tracks.pdf",
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
  .q <- vl_bw_average_track(bed= .SD, 
                            set_IDs = recovery,
                            tracks= track,
                            upstream = 2500,
                            downstream = 10000,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            legend= F,
                            col = Cc)
  
  # Legend
  leg <- unique(.q[, .(set_IDs, col)])
  leg[, legend("topright",
               legend = set_IDs,
               fill= col,
               bty= "n")]
  title(main= paste(variable, cdition))
  
  # Boxplot
  vl_boxplot(value~recovery,
             col= Cc,
             ylab= "Enrichment", 
             compute_pval= list(c(1,2)),
             tilt.names= T)
  print(".")
}, .(variable, cdition, track)]
dev.off()

#-----------------------------------------#
# Plot 2
#-----------------------------------------#
Cc <- c("chartreuse3", "brown2", "chocolate1", "darkorchid2")
sub <- dat[cdition %in% c("PH18", "PH29", "PHD9", "PHD11")]
sub[, cdition:= factor(cdition, c("PH18", "PH29", "PHD9", "PHD11"))]
setorderv(sub, "cdition")

pdf("pdf/recovery_cditions_average_tracks.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
sub[, {
  # Plot average track for each class
  .q <- vl_bw_average_track(bed= .SD,
                            tracks= unique(track),
                            upstream = 2500,
                            downstream = 10000,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            col = Cc,
                            legend= F)
  # Legend
  leg <- unique(.q[, .(name, col)])
  leg[, name:= tstrsplit(name, "_", keep= 2)]
  leg[, legend("topright",
               legend = name,
               fill= col,
               bty= "n",
               cex= 0.7)]
  
  title(main= paste(recovery, variable))
  
  # Boxplot
  vl_boxplot(value~cdition,
             col= adjustcolor(Cc, 0.5),
             ylab= "Enrichment", 
             compute_pval= list(c(1,2), c(1,3), c(1,4)),
             tilt.names= T)
  print(".")
}, .(recovery, variable)]
dev.off()

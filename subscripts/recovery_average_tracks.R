setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(rtracklayer)
require(GenomicRanges)
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")[, .(seqnames, start, end, strand, recovery)]
# Select recovery genes whose TSS is bound by PRC1
dat <- vl_resizeBed(dat[!is.na(recovery)], "start", 0, 0)
PRC1 <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
dat <- vl_intersectBed(dat, PRC1)
# bw files
tracks <- data.table(file= c(list.files("db/bw/cutnrun/", "merge", full.names = T),
                             "db/bw/ATAC/ATAC_merged.bw",
                             "db/bw/cutnrun_EcR/EcR_-6hAPF_merge.bw",
                             "db/bw/cutnrun_EcR/EcR_+6hAPF_merge.bw",
                             "db/bw/SA_2020/PC_ED_merge.bw",
                             "db/bw/SA_2020/PSC_ED_merge.bw",
                             "db/bw/SA_2020/PH_ED_merge.bw",
                             "db/bw/SA_2020/SUZ12_ED_merge.bw",
                             "db/bw/SA_2020/PHO_CNSID_merge.bw",
                             "db/bw/SA_2020/COMBGAP_CNSID_RP.bw",
                             "db/bw/SA_2020/SPPS_CNSID_BL_merge.bw",
                             "db/bw/SA_2020/GAF_WD_OH_merge.bw"))
tracks[, c("variable", "cdition"):= tstrsplit(basename(file), "_", keep= 1:2)]
tracks[variable=="EcR", c("variable", "cdition"):= .(paste0(variable, cdition), "ED")]
tracks[variable=="ATAC", cdition:= "ED"]
tracks[grepl("SA_2020", file), variable:= paste0(variable,"_SciAdv")]
# Gene body
tracks[grepl("^H", variable) & variable!="H3K27Ac", c("upstream", "downstream"):= .(1000, 5000)]
# Gene promoter
tracks[is.na(upstream), c("upstream", "downstream"):= .(750, 500)]
# Resize regions
dat <- tracks[, {
  vl_resizeBed(dat, 
               center = "start", 
               upstream = upstream, 
               downstream = downstream, 
               genome = "dm6")
}, (tracks)]
# Order before plotting
dat[, recovery:= factor(recovery, c("Recovery", "noRecovery"))]
dat[, cdition:= factor(cdition, c("PH18", "PH29", "PHD9", "PHD11", "ED", "CNSID", "WD"))]
setorderv(dat, c("recovery", "variable", "cdition"))

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
    tcl= -0.2,
    cex.lab= 0.9)
dat[!(cdition %in% c("PH29", "PHD9", "PHD11")), {
  # Plot average track for each class
  vl_bw_average_track(bed= .SD, 
                      set_IDs = recovery,
                      tracks= file,
                      upstream = 2500,
                      downstream = 10000,
                      stranded = T,
                      center_label = "TSS",
                      legend = F,
                      col = Cc)
  # Add quantification limits
  abline(v= c(-upstream, downstream), lty= 2)
  title(main= paste(variable, cdition))
  legend("topright", 
         legend= c("Recovery", "noRecovery"), 
         fill= Cc,
         bty= "n")
  # Boxplot
  ylab <- if(upstream==(1000))
    "Enrichment (TSS -1kb/+5kb)" else
      "Enrichment (TSS -750/+500bp)"
  value <- vl_bw_coverage(.SD, file)
  vl_boxplot(value~recovery,
             col= Cc,
             ylab= ylab, 
             compute_pval= list(c(1,2)),
             tilt.names= T)
  print(".")
}, .(variable, cdition, file, upstream, downstream)]
dev.off()

#-----------------------------------------#
# Plot 2
#-----------------------------------------#
Cc <- c("chartreuse3", "brown2", "chocolate1", "darkorchid2")

pdf("pdf/recovery_cditions_average_tracks.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
dat[cdition %in% c("PH18", "PH29", "PHD9", "PHD11"), {
  # Plot average track for each class
  vl_bw_average_track(bed= .SD,
                      tracks= unique(file),
                      upstream = 2500,
                      downstream = 10000,
                      stranded = T,
                      center_label = "TSS",
                      col = Cc,
                      legend= F)
  # Legend
  abline(v= c(-upstream, downstream), lty= 2)
  title(main= paste(recovery, variable))
  legend("topright", 
         legend= c("PH18", "PH29", "PHD9", "PHD11"), 
         fill= Cc,
         bty= "n")
  
  # Boxplot
  ylab <- if(upstream==(1000))
    "Enrichment (TSS -1kb/+5kb)" else
      "Enrichment (TSS -750/+500bp)"
  value <- lapply(unique(file), function(x) vl_bw_coverage(.SD, x))
  vl_boxplot(value,
             col= adjustcolor(Cc, 0.5),
             compute_pval= list(c(1,2), c(1,3), c(1,4)),
             tilt.names= T,
             names= levels(droplevels(cdition)),
             ylab= ylab)
  print(".")
}, .(recovery, variable, upstream, downstream)]
dev.off()

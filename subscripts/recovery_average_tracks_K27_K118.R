# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")[, .(seqnames, start, end, strand, recovery)]
dat <- dat[!is.na(recovery)]
# bw files
tracks <- data.table(file= c("db/bw/cutnrun/H3K27me3_PH18_merge.bw",
                             "db/bw/cutnrun/H2AK118Ub_PH18_merge.bw"))
tracks[, file:= factor(file, file)]
tracks[, variable:= gsub("_merge.bw$", "", basename(as.character(file)))]
# Gene body/promoter
tracks[, c("upstream", "downstream"):= .(1000, 5000)]
# Resize regions for quantification
dat <- tracks[, {
  vl_resizeBed(dat, 
               center = "start", 
               upstream = upstream, 
               # genome = "dm6",
               downstream = downstream) 
}, (tracks)]
# Order before plotting
dat[, recovery:= factor(switch(recovery, 
                               "noRecovery"= "Irreversible",
                               "Recovery"= "Reversible"), 
                        c("Irreversible", "Reversible")), recovery]
setorderv(dat, c("file", "recovery"))

#-----------------------------------------#
# Plot 1
#-----------------------------------------#
Cc <- c("rosybrown1", "palegreen3")
# Cc <- c("chartreuse3", "brown2", "chocolate1", "darkorchid2")
# Cc <- c("darkorchid2", "chartreuse3")

pdf("pdf/Extended_data_5_average_tracks.pdf",
    height = 6.5, 
    width = 35)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.4), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.35,0),
    tcl= -0.2,
    cex= 1)
dat[, {
  # Plot average track for each class
  vl_bw_average_track(bed= .SD, 
                      set_IDs = recovery,
                      tracks= as.character(file),
                      upstream = 2500,
                      downstream = 10000,
                      stranded = T,
                      center_label = "TSS",
                      legend= F,
                      col = Cc, 
                      ylab= paste(tstrsplit(variable, "_", keep= 1), "average signal"))
  # Add quantification limits
  abline(v= c(-upstream, downstream), lty= 2)
  title(main= "No ph-KD")
  legend("topright", 
         legend= c(paste0("Irreversible (n= ", sum(recovery=="Irreversible"), ")"), 
                   paste0("Reversible (n= ", sum(recovery=="Reversible"), ")")), 
         fill= adjustcolor(Cc, 0.6),
         bty= "n")
  # Boxplot
  ylab <- if(upstream==(1000))
    "Gene body enrich. ( -1kb/+5kb)" else
      "TSS enrich. (+/-500bp)"
  value <- vl_bw_coverage(.SD, as.character(file))
  vl_boxplot(value~recovery,
             col= adjustcolor(Cc, 0.6),
             ylab= ylab, 
             compute_pval= list(c(1,2)),
             tilt.names= T,
             srt= 30)
  print(".")
}, .(variable, file, upstream, downstream)]
dev.off()
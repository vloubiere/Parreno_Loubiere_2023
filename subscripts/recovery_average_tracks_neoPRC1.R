setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(rtracklayer)
require(GenomicRanges)
require(vlfunctions)

# Import data
SA2020 <- get(load("external_data/SA2020_cl.list"))
PH <- data.table(SA2020$ED_summits, SA2020$clusters_simp)
setnames(PH, c("seqnames", "start", "end", "name", "enrich", "recovery"))
neoPRC1 <- PH[recovery %in% 1:2]
neoPRC1[, recovery:= as.character(c("aTSS A+B", "enhancer A+B")[recovery])]
# Select PH peaks overlapping recovery/no recovery TSSs
recov <- fread("Rdata/final_gene_features_table.txt")
recov <- recov[!is.na(recovery), .(seqnames, start, end, strand, recovery)]
recov <- vl_resizeBed(recov, "start", 750, 500)
PH[, recovery:= as.character(recovery)]
PH[recov, recovery:= i.recovery, on= c("seqnames", "start<=end", "end>=start")]
recov <- PH[recovery %in% c("noRecovery", "Recovery")]
# object
dat <- rbind(recov,
             neoPRC1,
             fill= TRUE)
dat <- unique(dat[, .(seqnames, start, end, recovery)])
dat[, recovery:= factor(recovery, c("Recovery", "noRecovery", "aTSS A+B",  "enhancer A+B"))]
# Quantify signal
dat[, value:= vl_bw_coverage(vl_resizeBed(.SD, "center", 250, 250), 
                             bw = "db/bw/SA_2020/PH_ED_merge.bw")]
setorderv(dat, "recovery")

#-----------------------------------------#
# Plot 1
#-----------------------------------------#
Cc <- c("palegreen3", "rosybrown1", "tomato", "goldenrod")

pdf("pdf/recovery_neoPRC1_average_tracks.pdf",
    height = 4, 
    width = 7)
layout(matrix(1:2, ncol=2, byrow = T), 
       widths= c(1,0.5))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2,
    cex.lab= 0.9)
# Plot average track for each class
vl_bw_average_track(bed= dat, 
                    set_IDs = dat$recovery,
                    tracks= "db/bw/SA_2020/PH_ED_merge.bw",
                    upstream = 1000,
                    downstream = 1000,
                    stranded = F,
                    center_label = "PH ChIP-Seq peak",
                    legend = F,
                    col = Cc)
# Add quantification limits
abline(v= c(-250, 250), lty= 2)
title(main= "PH ChIP-Seq WT ED (SA 2020)")
legend("topright", 
       legend= levels(dat$recovery), 
       fill= Cc,
       bty= "n",
       cex= 0.9)
# Boxplot
vl_boxplot(value~recovery,
           dat,
           col= Cc,
           ylab= "Enrichment", 
           compute_pval= list(c(1,2), c(2,3), c(3,4)),
           tilt.names= T,
           cex= 0.9)
dev.off()
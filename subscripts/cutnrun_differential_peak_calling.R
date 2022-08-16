setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# Import peaks and K27 nonoverlapping peaks
peaks <- fread("Rdata/processed_metadata_CUTNRUN.txt")
peaks <- peaks[ChIP %in% c("H3K27me3", "H3K27Ac"), vl_importBed(peaks_merged), .(cdition, ChIP, peaks_merged)]
peaks <- peaks[V9>2 & V7>2]
nonoverlapping <- vl_importBed("db/narrowpeaks/cutnrun/K27_segmentation_merged.bed")

# upset plot
K27 <- split(peaks, list(peaks$ChIP, peaks$cdition))
K27 <- lapply(K27, function(x) {
  nonoverlapping[x[nonoverlapping, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N>0, paste0(seqnames, ":", start, "-", end)]
})

pdf("pdf/cutnrun/upset_plot_overlap_peaks_between_conditions.pdf", width = 10)
vl_upset_plot(K27,
              intersection_cutoff = 20)
dev.off()

# Alluvial plot
setkeyv(peaks, c("seqnames", "start", "end"))
setkeyv(nonoverlapping, c("seqnames", "start", "end"))
peaks[, {
  .c <- foverlaps(.SD, nonoverlapping)
  .c <- .c[, .(ChIP= ChIP[which.max(min(c(end, i.end))-max(c(start, i.start)))]), .(seqnames, start, end)]
  nonoverlapping[.c, (cdition):= i.ChIP, on= c("seqnames", "start", "end")]
  print("")
}, cdition]
cols <- grep("^PH", names(nonoverlapping))
nonoverlapping[, (cols):= lapply(.SD, function(x) {x[is.na(x)] <- "unmarked"; return(x)}), .SDcols= cols]
nonoverlapping <- nonoverlapping[!(PH18=="unmarked" 
                                   & PHD11=="unmarked"
                                   & PHD9=="unmarked"
                                   & PH29=="unmarked")]

pdf("pdf/cutnrun/alluvial_plot_nonOverlapping_K27_regions.pdf", width = 10)
vl_alluvial_plot(nonoverlapping[, .(PH18, PHD11, PHD9, PH29)])
dev.off()

setcolorder(nonoverlapping, c("PH18", "PHD11", "PHD9", "PH29"))
nonoverlapping[, c("group", "size"):= .(.GRP, .N), .(PH18, PHD11, PHD9, PH29)]
setorderv(nonoverlapping, c("PH18", "PHD11", "PHD9", "PH29"))

mat <- nonoverlapping[, .(coor= paste0(seqnames, ":", start, "-", end), PH18, PHD11, PHD9, PH29)]
cols <- grep("PH", names(mat))
mat[, (cols):= lapply(.SD, function(x) sapply(x, function(y) switch(y, "H3K27Ac"= 1, "H3K27me3"= -1, "unmarked"= 0))), .SDcols= cols]
vl_heatmap(as.matrix(mat[PH18!=0], 1), cluster_rows = F, cluster_cols = F)

vl_heatmap(as.matrix(mat[PH18==0], 1), cluster_rows = F, cluster_cols = F)



# nonoverlapping[PH18=="H3K27Ac" & PHD9=="unmarked" & PHD11=="unmarked" & PH29=="unmarked", K27Ac_class:= "ac_lost_all"]
# nonoverlapping[PH18=="H3K27Ac" & PHD9=="H3K27Ac" & PHD11=="H3K27Ac" & PH29=="unmarked", K27Ac_class:= "ac_lost_PH29"]
# nonoverlapping[PH18=="H3K27Ac" & PHD9=="H3K27Ac" & PHD11=="unmarked" & PH29=="H3K27Ac", K27Ac_class:= "ac_lost_PHD11"]
# nonoverlapping[PH18=="H3K27Ac" & PHD9=="unmarked" & PHD11=="H3K27Ac" & PH29=="H3K27Ac", K27Ac_class:= "ac_lost_PHD9"]
# nonoverlapping[PH18=="unmarked" & PHD9=="H3K27Ac" & PHD11=="H3K27Ac" & PH29=="H3K27Ac", K27Ac_class:= "ac_gain_all"]
# nonoverlapping[PH18=="unmarked" & PHD9=="unmarked" & PHD11=="unmarked" & PH29=="H3K27Ac", K27Ac_class:= "ac_gain_PH29"]
# nonoverlapping[PH18=="unmarked" & PHD9=="unmarked" & PHD11=="H3K27Ac" & PH29=="unmarked", K27Ac_class:= "ac_gain_PHD11"]
# nonoverlapping[PH18=="unmarked" & PHD9=="H3K27Ac" & PHD11=="unmarked" & PH29=="unmarked", K27Ac_class:= "K27_gain_PHD9"]
# 
# nonoverlapping[PH18=="H3K27me3" & PHD9=="unmarked" & PHD11=="unmarked" & PH29=="unmarked", K27me3_class:= "me3_lost_all"]
# nonoverlapping[PH18=="H3K27me3" & PHD9=="H3K27me3" & PHD11=="H3K27me3" & PH29=="unmarked", K27me3_class:= "me3_lost_PH29"]
# nonoverlapping[PH18=="H3K27me3" & PHD9=="H3K27me3" & PHD11=="unmarked" & PH29=="H3K27me3", K27me3_class:= "me3_lost_PHD11"]
# nonoverlapping[PH18=="H3K27me3" & PHD9=="unmarked" & PHD11=="H3K27me3" & PH29=="H3K27me3", K27me3_class:= "me3_lost_PHD9"]
# nonoverlapping[PH18=="unmarked" & PHD9=="H3K27me3" & PHD11=="H3K27me3" & PH29=="H3K27me3", K27me3_class:= "me3_gain_all"]
# nonoverlapping[PH18=="unmarked" & PHD9=="unmarked" & PHD11=="unmarked" & PH29=="H3K27me3", K27me3_class:= "me3_gain_PH29"]
# nonoverlapping[PH18=="unmarked" & PHD9=="unmarked" & PHD11=="H3K27me3" & PH29=="unmarked", K27me3_class:= "me3_gain_PHD11"]
# nonoverlapping[PH18=="unmarked" & PHD9=="H3K27me3" & PHD11=="unmarked" & PH29=="unmarked", K27me3_class:= "me3_gain_PHD9"]





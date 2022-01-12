setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
options(scipen = 999)
require(kohonen)
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm6)

reads <- fread("Rdata/metadata_cutnrun_final.txt")
reads[, shuffled_bed:= paste0("db/bed/cutnrun/merge/", ChIP, "_", cdition, "_merge_shuffled.bed"), .(ChIP, cdition)]
reads[, peaks_narrowpeaks:= paste0("db/narrowpeaks/", ChIP, "_", cdition, "_peaks.narrowPeak"), .(ChIP, cdition)]
reads[, peaks_merged:= paste0("db/narrowpeaks/", ChIP, "_", cdition, "_merged_peaks.narrowPeak"), .(ChIP, cdition)]

# Peaks calling
reads[, {
  bs <- switch(ChIP, "H3K27Ac"= 100, "H3K27me3"= 1000)
  min_width <- switch(ChIP, "H3K27Ac"= 250, "H3K27me3"= 5000)
  min_dist <- switch(ChIP, "H3K27Ac"= 201, "H3K27me3"= 2501)
  min_qval <- switch(ChIP, "H3K27Ac"= 3, "H3K27me3"= 10)
  min_OR <- switch(ChIP, "H3K27Ac"= 1, "H3K27me3"= 2)
  # Create shuffled bed if not already done
  if(!file.exists(shuffled_bed))
    vl_exportBed(vl_shuffleBed(ChIP_bed_merge), shuffled_bed)
  # Peak calling
  if(!file.exists(peaks_narrowpeaks))
  {
    peaks <- vl_peakCalling(ChIP_bed_merge,
                            shuffled_bed,
                            gaussian_blur = T,
                            BSgenome = BSgenome.Dmelanogaster.UCSC.dm6, 
                            bins_width = bs, 
                            bins_OR_cutoff = 1)
    fwrite(peaks, 
           peaks_narrowpeaks, 
           col.names = F, 
           sep= "\t")
  }
  # Merged peaks
  if(!file.exists(peaks_merged))
  {
    if(!exists(peaks_narrowpeaks))
      peaks <- vl_importBed(peaks_narrowpeaks)
    .m <- vl_collapseBed(peaks, mingap = min_dist)[end-start>min_width]
    merged <- vl_enrichBed(.m, 
                           ChIP_bed_merge, 
                           shuffled_bed)
    fwrite(merged[qValue>min_qval & signalValue>=min_OR], 
           peaks_merged,
           col.names = F, 
           sep= "\t")
  }
  print("done")
}, .(ChIP, cdition, ChIP_bed_merge, shuffled_bed, peaks_narrowpeaks, peaks_merged)]

# Segment genome into homogenous regions
if(!file.exists("Rdata/K27_regions_segmentation.bed"))
{
  merged_peaks <- reads[, vl_importBed(peaks_merged), .(ChIP, cdition, peaks_merged)]
  merged_peaks <- vl_collapseBed(merged_peaks, return_idx_only = T)
  setkeyv(merged_peaks, c("seqnames", "start"))
  K27 <- merged_peaks[, {
    ranges <- sort(unique(c(start, end[-.N]+1, end[.N])))
    .(start= ranges[-length(ranges)], end= ranges[-1]-1)
  }, .(seqnames, idx)]
  K27 <- unique(K27)[end-start>300]
  K27 <- cbind(K27,
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_shuffled.bed")[, .(K27Ac_WT= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PHD11_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed")[, .(K27Ac_PHD11= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PHD9_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed")[, .(K27Ac_PHD9= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PH29_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed")[, .(K27Ac_PH29= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_shuffled.bed")[, .(K27me3_WT= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PHD11_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed")[, .(K27me3_PHD11= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PHD9_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed")[, .(K27me3_PHD9= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PH29_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed")[, .(K27me3_PH29= signalValue)])
  fwrite(na.omit(K27),
         "Rdata/K27_regions_segmentation.bed",
         sep= "\t", 
         col.names = T)
}

#-------------------#
# clustering
#-------------------#
# Import
dat <- fread("Rdata/K27_regions_segmentation.bed")
cols <- grep("^K27", names(dat), value = T)
dat[, (cols):= lapply(.SD, function(x) {
  lim <- quantile(x, c(0.001, 0.999))
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  log2(x)
}), .SDcols= cols]
# blacklisted <- vl_importBed("/mnt/d/_R_data/genomes/dm6/dm6-blacklist.v2.bed")
# Clip extremes
grid <- somgrid(xdim = 2, 
                ydim = 3, 
                topo = "hexagonal", 
                toroidal = T)
som <- supersom(list(as.matrix(dat[, .(K27Ac_WT, K27me3_WT)]),
                     as.matrix(dat[, .(K27Ac_PHD11, K27Ac_PHD9, K27Ac_PH29, K27me3_PHD11, K27me3_PHD9, K27me3_PH29)])), 
                grid = grid, 
                user.weights = c(10, 4))
              
dat[, cl:= som$unit.classif]
setkeyv(dat, "cl")

par(mfrow= c(1,2))
setcolorder(dat, 
            c("cl", "K27Ac_WT", "K27me3_WT"))
vl_heatmap(dat[, cl:K27me3_WT], 
           cluster_rows= F,
           cluster_cols= F, 
           breaks= c(-2,-0.25,0.25,2), 
           col= c("cornflowerblue", "white", "white", "tomato"))
setcolorder(dat, 
            c("cl", "K27Ac_PHD11", "K27Ac_PHD9", "K27Ac_PH29", "K27me3_PHD11", "K27me3_PHD9", "K27me3_PH29"))
vl_heatmap(dat[, cl:K27me3_PH29], 
           cluster_rows= F,
           cluster_cols= F, 
           auto_margins= F, 
           breaks= c(-2,-0.25,0.25,2), 
           col= c("cornflowerblue", "white", "white", "tomato"))





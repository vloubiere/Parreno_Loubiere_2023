setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
options(scipen = 999)
require(kohonen)
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm6)

reads <- fread("Rdata/processed_metadata_CUTNRUN.txt")
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

#--------------------------------------------#
# Segment genome into homogenous regions
#--------------------------------------------#
# Import peaks
merged_peaks <- reads[, vl_importBed(peaks_merged), .(ChIP, cdition, peaks_merged)]
merged_peaks <- vl_collapseBed(merged_peaks, return_idx_only = T)
setkeyv(merged_peaks, c("seqnames", "start"))
# Merge homogneous blocks
K27 <- merged_peaks[, {
  ranges <- sort(unique(c(start, end[-.N]+1, end[.N])))
  .(start= ranges[-length(ranges)], end= ranges[-1]-1)
}, .(seqnames, idx)]
K27 <- unique(K27)[end-start>300]
# Compute enrichment over 18C control
ChIP <- reads[, cbind(K27,
                      counts= vl_covBed(K27, ChIP_bed_merge),
                      total_counts= fread(cmd= paste0("wc -l ", ChIP_bed_merge))$V1), .(ChIP_bed_merge, ChIP, cdition)]
shuffle <- reads[cdition=="PH18", cbind(K27,
                                        counts= vl_covBed(K27, shuffled_bed),
                                        total_counts= fread(cmd= paste0("wc -l ", shuffled_bed))$V1), .(shuffled_bed, ChIP, cdition)]
# WT levels (over shuffled)
res <- rbind(merge(ChIP[cdition=="PH18", .(ChIP, cdition, seqnames, start, end, counts, total_counts)],
                   shuffle[cdition=="PH18", .(ChIP, seqnames, start, end, counts, total_counts)],
                   by= c("ChIP", "seqnames", "start", "end"),
                   suffixes= c("_ChIP", "_control")),
             # Mutant levels over PH18
             merge(ChIP[cdition!="PH18", .(ChIP, cdition, seqnames, start, end, counts, total_counts)],
                   ChIP[cdition=="PH18", .(ChIP, seqnames, start, end, counts, total_counts)],
                   by= c("ChIP", "seqnames", "start", "end"),
                   suffixes= c("_ChIP", "_control")))
check <- res[,counts_ChIP>0 & counts_control>0]
res[(check), c("OR", "pval"):= {
  mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
  fisher.test(mat)[c("estimate", "p.value")]
}, .(counts_ChIP, counts_control, total_counts_ChIP, total_counts_control)]
res[, padj:= p.adjust(pval, method= "fdr")]
saveRDS(res, "Rdata/K27_regions_segmentation.rds")
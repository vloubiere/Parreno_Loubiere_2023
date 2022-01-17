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

# Import peaks
merged_peaks <- reads[, vl_importBed(unique(peaks_merged)), ChIP]
# Compute K27me3 counts
K27me3 <- vl_collapseBed(merged_peaks[ChIP=="H3K27me3"])
K27me3 <- reads[ChIP=="H3K27me3", cbind(K27me3,
                                        counts= vl_covBed(K27me3, ChIP_bed_merge),
                                        total_counts= fread(cmd= paste0("wc -l ", ChIP_bed_merge))$V1), .(ChIP_bed_merge, ChIP, cdition)]
K27me3 <- merge(K27me3[, .(ChIP, cdition, seqnames, start, end, counts, total_counts)],
                K27me3[cdition=="PH18", .(ChIP, seqnames, start, end, counts, total_counts)],
                by= c("ChIP", "seqnames", "start", "end"),
                suffixes= c("_ChIP", "_control"))
# Compute K27Ac counts
K27Ac <- vl_collapseBed(merged_peaks[ChIP=="H3K27Ac"])
K27Ac <- reads[ChIP=="H3K27Ac", cbind(K27Ac,
                                      counts= vl_covBed(K27Ac, ChIP_bed_merge),
                                      total_counts= fread(cmd= paste0("wc -l ", ChIP_bed_merge))$V1), .(ChIP_bed_merge, ChIP, cdition)] 
K27Ac <- merge(K27Ac[, .(ChIP, cdition, seqnames, start, end, counts, total_counts)],
               K27Ac[cdition=="PH18", .(ChIP, seqnames, start, end, counts, total_counts)],
               by= c("ChIP", "seqnames", "start", "end"),
               suffixes= c("_ChIP", "_control"))

# FC over PH18
res <- rbind(K27me3,
             K27Ac)
check <- res[,counts_ChIP>0 & counts_control>0]
res[(check), c("OR", "pval"):= {
  mat <- matrix(unlist(.BY), nrow= 2, byrow = T)
  fisher.test(mat)[c("estimate", "p.value")]
}, .(counts_ChIP, counts_control, total_counts_ChIP, total_counts_control)]

res[, padj:= p.adjust(pval, method= "fdr")]


res[, ChIP:= factor(ChIP, levels= c("H3K27me3", "H3K27Ac"))]
res[, cdition:= factor(cdition, levels= c("PH18", "PHD11", "PHD9", "PH29"))]

par(mfrow=c(2,4))
res[, plot(log2(OR), 
           -log10(padj),
           main= paste(ChIP, cdition),
           xlim= switch(as.character(ChIP), 
                        "H3K27Ac"= c(-3,3),
                        "H3K27me3"= c(-2,2)),
           ylim= c(0, switch(as.character(ChIP), 
                             "H3K27Ac"= 200,
                             "H3K27me3"= 350))), keyby= .(ChIP, cdition)]

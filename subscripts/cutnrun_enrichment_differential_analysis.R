setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
options(scipen = 999)
require(kohonen)
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm6)

reads <- fread("Rdata/processed_metadata_CUTNRUN.txt")

# K27me3 peaks
merged_peaks <- reads[, vl_importBed(unique(peaks_merged)), ChIP]
# Compute K27me3 counts
K27me3 <- vl_collapseBed(merged_peaks[ChIP=="H3K27me3"])
K27me3 <- reads[ChIP=="H3K27me3", cbind(K27me3,
                                        counts= vl_covBed(K27me3, ChIP_bed)), .(ChIP_bed, ChIP, cdition, rep)]
K27me3[, coor:= paste0(seqnames, ":",  start, "-", end)]
DF <- dcast(K27me3, coor~cdition+rep, value.var = "counts")
DF <- data.frame(DF[, -1], row.names = DF$coor)
sampleTable <- data.table(colnames(DF))
sampleTable[, c("cdition", "rep"):= tstrsplit(V1, "_")]
sampleTable <- data.frame(sampleTable[,-1], 
                          row.names = sampleTable[,V1])

dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                      colData= sampleTable,
                                      design= ~rep+cdition)















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
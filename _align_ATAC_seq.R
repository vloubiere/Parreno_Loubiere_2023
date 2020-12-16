setwd("/_R_data/projects/epigenetic_cancer/")
require(Rsubread)
require(DESeq2)
require(data.table)
require(BSgenome.Dmelanogaster.UCSC.dm6)

#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
dat <- data.table(file= list.files("../../public_data/dm6/fq/ATAC_JJ_2018/", full.names = T))
dat[, bam:= paste0("../../public_data/dm6/bam/ATAC_JJ_2018/", gsub(".fastq.gz", ".bam", basename(file)))]
dat[, {
  if(!file.exists(bam)){
    stats <- capture.output(align(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", readfile1= file, nthreads = 8, unique = T, output_file= bam))
    # align(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", readfile1= file, nthreads = 8, unique = T, output_file= bam)
    writeLines(stats, con = gsub(".bam$", "_stats.txt", bam))
  }
  print(paste(bam, "DONE!"))
}, .(file, bam)]

#----------------------------------------------------------#
# Compute counts
#----------------------------------------------------------#
if(!file.exists("../../public_data/dm6/counts/ATAC_gw_bin_counts.rds")){
  chr <- as.data.table(seqinfo(BSgenome.Dmelanogaster.UCSC.dm6), keep.rownames= "Chr")
  bins <- chr[Chr %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"), .(Start= seq(1, seqlengths, 10)), .(Chr, seqlengths)]
  bins[, End:= ifelse(Start+9>seqlengths, seqlengths, Start+9)]
  bins[, Strand:= "+"]
  bins[, GeneID:= paste0(Chr, ":", Start, "-", End)]
  bins <- bins[, .(GeneID, Chr, Start, End, Strand)]
  counts <- featureCounts(dat$bam, annot.ext= as.data.frame(bins), isPairedEnd = F, nthreads = 8, allowMultiOverlap = T)
  saveRDS(counts, "../../public_data/dm6/counts/ATAC_gw_bin_counts.rds")
}

#----------------------------------------------------------#
# peak calling
#----------------------------------------------------------#
if(!file.exists("../../public_data/dm6/peaks/merged_ATAC_seq_peaks.txt")){
  counts <- as.data.table(readRDS("../../public_data/dm6/counts/ATAC_gw_bin_counts.rds")$counts, keep.rownames= T)
  # COmpute enrichment and filter peaks that are consistently found across strains
  .m <- melt(counts, id.vars = "rn")
  .m[, norm:= log2(value+1)-log2(mean(value)+1), variable]
  peaks <- .m[norm>3, .(N= .N), rn]
  peaks <- GRanges(peaks[N>12, rn]) #enriched in at leat 12/16 conditions
  # Merge peaks and quantify signal
  peaks <- as.data.table(reduce(peaks, min.gapwidth= 101))
  peaks <- peaks[, .(GeneID= paste0(seqnames, ":", start, "-", end), Chr= seqnames, Start= start, End= end, Strand= "+")]
  enrich <- as.data.table(featureCounts(list.files("D:/_R_data/public_data/dm6/bam/ATAC_JJ_2018/", ".bam$", full.names = T), 
                                        annot.ext= as.data.frame(peaks), isPairedEnd = F, nthreads = 8, allowMultiOverlap = T)$counts)
  peaks[, log2_counts:= log2(rowSums(enrich)+1)]
  fwrite(peaks, "../../public_data/dm6/peaks/merged_ATAC_seq_peaks.txt", col.names = T, row.names = F, sep= "\t", quote= F)
  export(peaks, "../../public_data/dm6/peaks/merged_ATAC_seq_peaks.bed")
}

#----------------------------------------------------------#
# bedgraph generation
#----------------------------------------------------------#
if(!file.exists("../../public_data/dm6/bg/ATAC_JJ_ED_merge.bedgraph")){
  if(!exists("counts")){
  counts <- as.data.table(readRDS("../../public_data/dm6/counts/ATAC_gw_bin_counts.rds")$counts, keep.rownames= T)
  }
  bdg <- GRanges(counts$rn)
  bdg$score <- counts[, rowSums(.SD), .SDcols= patterns("^SRR")]
  export.bedGraph(bdg, "../../public_data/dm6/bg/ATAC_JJ_ED_merge.bedgraph")
}



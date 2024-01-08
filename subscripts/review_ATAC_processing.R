# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(parallel)
require(readxl)

meta <- fread("Rdata/metadata_ATAC_review.txt")
meta[, fq1:= list.files("db/fq/ATAC", paste0(fq1, "$"), recursive = T, full.names = T), fq1]
meta[, fq2:= list.files("db/fq/ATAC", paste0(fq2, "$"), recursive = T, full.names = T), fq2]
meta[, fq1_trim:= gsub(".fq.gz$", "_val_1.fq.gz", fq1)]
meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]
meta[, bam_file:= paste0("/scratch/stark/vloubiere/bam/ATAC/", gsub("_1.fq.gz", ".bam", basename(fq1))), fq1]
meta[, bam:= basename(bam_file)]
meta[, peaks_reps:= paste0("db/peaks/ATAC/ATAC_", cdition, "_rep", rep, "_peaks.narrowPeak"), .(cdition, rep)]
meta[, peaks_merge:= paste0("db/peaks/ATAC/ATAC_", cdition, "_merge_peaks.narrowPeak"), cdition]
meta[, peaks_conf:= paste0("db/peaks/ATAC/ATAC_", cdition, "_conf_peaks.narrowPeak"), cdition]
meta[, bdg_reps:= paste0("db/peaks/ATAC/ATAC_", cdition, "_rep", rep, "_treat_pileup.bdg"), .(cdition, rep)]
meta[, bdg_merge:= paste0("db/peaks/ATAC/ATAC_", cdition, "_merge_treat_pileup.bdg"), cdition]
meta[, bw_reps:= paste0("db/bw/ATAC/ATAC_", cdition, "_rep", rep, ".bw"), .(cdition, rep)]
meta[, bw_merge:= paste0("db/bw/ATAC/ATAC_", cdition, "_merge.bw"), cdition]
meta[, saf_file:= paste0("db/saf/ATAC_peaks_", DESeq2, ".saf"), DESeq2]
meta[, counts_file:= paste0("db/counts/ATAC/ATAC_peaks_", DESeq2, "_ATAC_counts.rds"), DESeq2]
meta[, dds_file:= paste0("db/dds/ATAC/ATAC_FC_", DESeq2, "_ATAC_counts.dds"), DESeq2]

# Save processed metadata
fwrite(meta,
       "Rdata/metadata_ATAC_review_processed.txt",
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)

# Retrieve fastqs and trim reads
if(!all(file.exists(meta$fq1_trim, meta$fq2_trim)))
{
  mcmapply(FUN = function(outdir, fq1, fq2, fq1_trim, fq2_trim){
    if(!all(file.exists(c(fq1_trim, fq2_trim))))
    {
      cmd <- "module load build-env/2020; module load trim_galore/0.6.0-foss-2018b-python-2.7.15; "
      cmd <- paste0(cmd, "trim_galore --paired --gzip -o ", outdir, "/ ", fq1, " ", fq2)
      system(cmd)
    }
    print("done")
  },
  outdir= dirname(meta$fq1),
  fq1= meta$fq1,
  fq2= meta$fq2,
  fq1_trim= meta$fq1_trim,
  fq2_trim= meta$fq2_trim,
  mc.preschedule = F,
  mc.cores = getDTthreads())
}

# Align ----
meta[, {
  if(!file.exists(bam_file))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam_file)
    cmd <- paste0("module load build-env/2020; module load bowtie2/2.3.5.1-foss-2018b; bowtie2 -p 6 -x /groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1_trim, " -2 ", fq2_trim)
    cmd <- paste0(cmd, " -S ", sam_file)
    # sam to bam
    cmd <- paste0(cmd, "; module load samtools/0.1.19-foss-2018b; samtools view -@ 5 -S -b -q 30 ", sam_file, " -o ", bam_file)
    cmd <- paste0(cmd, "; rm ", sam_file)
    system(cmd)
    # vl_bsub(cmd, cores = 6, m= 16, t= '08:00:00') # Much faster if needed
  }
  print("DONE")
}, .(fq1_trim, fq2_trim, bam_file)]

# Peak calling replicates ----
meta[, {
  if(!file.exists(peaks_reps))
  {
    sorted <- gsub(".bdg$", "_s.bdg", bdg_reps)
    cmd <- paste0("module load build-env/2020; module load macs2/2.2.5-foss-2018b-python-3.6.6; macs2 callpeak --keep-dup 1 -g dm -f BAM --SPMR -B",
                  " --outdir db/peaks/ATAC/ ",
                  " -t ", bam_file,
                  paste0(" -n ", gsub("_peaks.narrowPeak", "", basename(peaks_reps))),
                  "; sort -k1,1 -k2,2n ", bdg_reps, " > ", sorted, " ;mv ", sorted, bdg_reps)
    system(cmd)
    # vl_bsub(cmd, cores = 4, m= 8, t= '08:00:00') # Much faster if needed
  }
  print("DONE")
}, .(peaks_reps, bdg_reps)]

# Peak calling merge ----
meta[, {
  if(!file.exists(peaks_merge))
  {
    sorted <- gsub(".bdg$", "_s.bdg", bdg_merge)
    cmd <- paste0("module load build-env/2020; module load macs2/2.2.5-foss-2018b-python-3.6.6; macs2 callpeak --keep-dup 1 -g dm -f BAM --SPMR -B",
                  " --outdir db/peaks/ATAC/ ",
                  " -t ", paste(bam_file, collapse = " "),
                  paste0(" -n ", gsub("_peaks.narrowPeak", "", basename(peaks_merge))),
                  "; sort -k1,1 -k2,2n ", bdg_merge, " > ", sorted, " ;mv ", sorted, bdg_merge)
    system(cmd)
    # vl_bsub(cmd, cores = 4, m= 8, t= '08:00:00') # Much faster if needed
  }
  print("DONE")
}, .(peaks_merge, bdg_merge)]

# Confident peaks ----
meta[, {
  if(!exists(peaks_conf))
  {
    conf <- vl_importBed(peaks_merge)[signalValue>2 & qValue>2]
    rep_peaks <- rbindlist(lapply(peaks_reps, vl_importBed, extraCols= "narrowPeak"), idcol = T)
    conf <- conf[rep_peaks[conf, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N==2]
    fwrite(conf[order(seqnames, start)],
           peaks_conf,
           col.names = F,
           sep= "\t",
           quote= F,
           na= NA)
  }
  print("DONE")
}, .(peaks_merge, peaks_conf)]

# bw files ----
meta[, {
  if(!file.exists(bw_reps))
  {
    cmd <- "/software/2020/software/kent_tools/20190507-linux.x86_64/bin/bedGraphToBigWig"
    cmd <- paste(cmd, bdg_reps, "/groups/stark/vloubiere/genomes/chromSizes/dm6_chromsizes.txt", bw_reps)
    system(cmd)
  }
  print("DONE")
}, .(bdg_reps, bw_reps)]
meta[, {
  if(!file.exists(bw_merge))
  {
    cmd <- "/software/2020/software/kent_tools/20190507-linux.x86_64/bin/bedGraphToBigWig"
    cmd <- paste(cmd, bdg_merge, "/groups/stark/vloubiere/genomes/chromSizes/dm6_chromsizes.txt", bw_merge)
    system(cmd)
  }
  print("DONE")
}, .(bdg_merge, bw_merge)]

# Compute features SAF files ----
meta[, {
  if(!file.exists(saf_file))
  {
    peaks <- .SD[, vl_importBed(peaks_conf), peaks_conf]
    merged <- vl_collapseBed(peaks, min.gap = 250)[, .(seqnames, start, end)]
    merged <- merged[, .(GeneID= paste0(seqnames, ":", start, "-", end),
                         Chr= seqnames,
                         Start= start,
                         End= end,
                         Strand= "*")]
    fwrite(merged,
           saf_file,
           sep= "\t")
  }
}, saf_file]

# Compute counts ----
meta[, {
  if(!file.exists(counts_file))
  {
    saf <- fread(saf_file)
    saf <- saf[Chr %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")]
    saf <- as.data.frame(saf)
    .c <- Rsubread::featureCounts(bam_file, # count reads
                                  annot.ext= saf,
                                  isGTFAnnotationFile = F,
                                  isPairedEnd = T,
                                  nthreads = 8)
    saveRDS(.c, counts_file)
  }
}, .(saf_file, counts_file)]

# DESeq2 analysis ----
meta[, {
  if(!file.exists(dds_file))
  {
    .c <- readRDS(counts_file)
    DF <- as.data.frame(.c[[1]])
    colnames(DF) <- paste0(cdition, "_", rep)
    DF <- DF[rowSums(DF[,-1])>100,]

    # Samples
    sampleTable <- data.table(sample= names(DF))
    sampleTable[, cdition:= gsub("(.*)_.*$", "\\1", sample)]
    sampleTable[, rep:= gsub(".*_(.*)$", "\\1", sample)]

    # DESeq2 dataset
    .dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                           colData= sampleTable,
                                           design= ~rep+cdition)
    # Libsize norm
    libsize <- .c[[4]][, -1]
    names(libsize) <- meta[names(libsize), paste0(cdition, "_", rep), on= "bam"]
    libsize <- apply(libsize, 2, sum)
    sizeFactors(.dds) <- libsize/min(libsize)
    # Run DESeq2 and save
    .dds <- DESeq2::DESeq(.dds)
    saveRDS(.dds, dds_file)
  }
}, .(dds_file, counts_file)]

# FC tables
FC1 <- data.table(V1= c("PH29", "PHD11", "PHD11"),
                  V2= c("PH18", "PH18", "PH29"), 
                  dds_file= "db/dds/ATAC/ATAC_FC_tumor_ATAC_counts.dds")
FC2 <- data.table(V1= c("GFP_PH", "Stat_PH", "Stat_W", "ZFH1_PH", "ZFH1_W"),
                  V2= "GFP_W", 
                  dds_file= "db/dds/ATAC/ATAC_FC_rescue_ATAC_counts.dds")
FC3 <- data.table(V1= c("Stat_PH", "ZFH1_PH"),
                  V2= "GFP_PH", 
                  dds_file= "db/dds/ATAC/ATAC_FC_rescue_ATAC_counts.dds")
FC <- rbind(FC1, FC2, FC3)
FC <- FC[V1!=V2]
FC[, output:= paste0("db/FC_tables/ATAC/ATAC_", V1, "_vs_", V2, ".txt")]
FC[, {
  if(!file.exists(output))
  {
    .dds <- readRDS(dds_file)
    res <- DESeq2::lfcShrink(.dds,
                             type= "ashr",
                             contrast= c("cdition",  V1,  V2))
    # res <- as.data.frame(DESeq2::results(.dds,
    #                                      contrast= c("cdition", V1, V2)))
    res <- as.data.frame(res)
    res <- as.data.table(res, keep.rownames = "ID")
    fwrite(res,
           output,
           col.names = T,
           sep= "\t",
           quote=F,
           na= NA)
  }
  print("DONE")
}, .(V1, V2, dds_file, output)]

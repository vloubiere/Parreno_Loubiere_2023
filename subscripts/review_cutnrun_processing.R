setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(data.table)
require(parallel)
require(Rsubread)
require(seqinr)
require(rtracklayer)

# METADATA -----
meta <- fread("Rdata/metadata_cutnrun_review.txt")
meta[, fq1:= list.files("db/fq/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = paste0(fq1, "$")), fq1]
meta[, fq2:= list.files("db/fq/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = paste0(fq2, "$")), fq2]
meta[, fq1_trim:= gsub(".fq.gz$", "_val_1.fq.gz", fq1)]
meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]
meta[, bam:= paste0("/scratch/stark/vloubiere/epicancer/", method, "/", ChIP, "_", cdition, "_", rep, ".bam")]
meta[meta[ChIP=="input"], input_bam:= i.bam, on= c("cdition", "rep", "method")]
meta[is.na(input_bam), input_bam:= "/scratch/stark/vloubiere/epicancer/input_PH18_1.bam"] # ZFH!/STAT92E had no Input so I used this one
meta[ChIP!="input", peaks_rep:= paste0("db/peaks/", method, "/", paste0(ChIP, "_", cdition, "_", rep))]
meta[, broad:= F] # Is it a ChIP/cutNrun for a broad mark
meta[ChIP!="input", peaks_merge:= paste0("db/peaks/", method, "/", paste0(ChIP, "_", cdition, "_merge"))]
meta[ChIP!="input", filtered_peaks:= paste0("db/peaks/", method, "/", ChIP, "_", cdition, 
                                            ifelse(broad, 
                                                   "_confident_peaks.broadPeak", 
                                                   "_confident_peaks.narrowPeak"))]
meta[, dist_cutoff:= 250]
meta[ChIP!="input", merged_file:= paste0("db/peaks/", method, "_merged_peaks/", ChIP, "_merged_peaks", 
                                         ifelse(broad, ".broadPeak", ".narrowPeak"))]
meta[ChIP!="input", bw_file:= paste0("db/bw/", method, "/", ChIP, "_", cdition, "_", rep, ".bw")]
meta[ChIP!="input", bw_merge:= paste0("db/bw/", method, "/", ChIP, "_", cdition, "_merge.bw")]
meta[ChIP!="input", peaks_saf:= paste0("db/saf/", ChIP, "_", "peaks.saf"), ChIP]
meta[ChIP!="input", peaks_counts:= paste0("db/counts/", method, "/",  ChIP, "_peaks_counts.txt")]

# meta <- meta[method=="ChIP"]

# Retrieve fastqs and trim reads ----
if(!all(file.exists(meta$fq1_trim, meta$fq2_trim)))
  mcmapply(FUN = function(outdir, fq1, fq2, fq1_trim, fq2_trim){
    if(!all(file.exists(c(fq1_trim, fq2_trim))))
      system(paste0("module load build-env/2020; module load trim_galore/0.6.0-foss-2018b-python-2.7.15; trim_galore --paired --gzip -o ",
                    outdir, "/ ", fq1, " ", fq2))
    print("done")
  },
  outdir= dirname(meta$fq1),
  fq1= meta$fq1,
  fq2= meta$fq2,
  fq1_trim= meta$fq1_trim,
  fq2_trim= meta$fq2_trim,
  mc.preschedule = F,
  mc.cores = getDTthreads())

# Alignment ----
meta[, {
  if(!file.exists(bam))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam)
    cmd <- paste0("module load build-env/2020; module load bowtie2/2.3.5.1-foss-2018b; bowtie2 -p 10 -x /groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1_trim, " -2 ", fq2_trim)
    cmd <- paste0(cmd, " -S ", sam_file)
    # sam to bam
    cmd <- paste0(cmd, "; module load samtools/0.1.19-foss-2018b; samtools view -@ 9 -S -b -q 30 ", sam_file, " -o ", bam)
    system(cmd)
    file.remove(sam_file)
  }
  print("DONE")
}, bam]

# Peak calling ----
MACS <- function(bam_ChIP,
                 bam_Input,
                 output,
                 broad,
                 SPMR)
{
  cmd <- paste0("module load build-env/2020; module load macs2/2.2.5-foss-2018b-python-3.6.6; macs2 callpeak -g dm --keep-dup 1 -f BAMPE",
                ifelse(SPMR, " -B --SPMR", ""),
                " --outdir ", paste0(dirname(output), "/"),
                " -t ", paste(bam_ChIP, collapse= " "), # Maybe add unique()?
                " -c ", paste(bam_Input, collapse= " "), # Maybe add unique()?
                " -n ", basename(output))
  if(broad)
    cmd <- paste0(cmd, " --broad")
  return(cmd)
}
## Separated replicates ----
meta[!is.na(peaks_rep), {
  if(!file.exists(paste0(peaks_rep, "_peaks.xls")))
  {
    cmd <- MACS(bam, 
                input_bam, 
                peaks_rep, 
                broad = broad,
                SPMR= T)
    system(cmd)
    print("DONE")
  }
}, .(ChIP, peaks_rep, broad)]
## Merge ----
meta[!is.na(peaks_merge), {
  if(!file.exists(paste0(peaks_merge, "_peaks.xls")))
  {
    cmd <- MACS(bam, 
                input_bam, 
                peaks_merge, 
                broad = broad,
                SPMR= F)
    system(cmd)
    print("DONE")
  }
}, .(ChIP, peaks_merge, broad)]
## Return peaks file ----
meta[!is.na(peaks_rep), peaks_rep:= paste0(peaks_rep, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]
meta[!is.na(peaks_merge), peaks_merge:= paste0(peaks_merge, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]

# Confident peaks ----
meta[!is.na(filtered_peaks), {
  if(!file.exists(filtered_peaks))
  {
    .c <- vl_importBed(peaks_merge)
    .c <- .c[signalValue>2 & qValue>2]
    if(any(rep==1))
      .c <- .c[vl_covBed(.c, unique(peaks_rep[rep==1]))>0]
    if(any(rep==2))
      .c <- .c[vl_covBed(.c, unique(peaks_rep[rep==2]))>0]
    .c$idx <- vl_collapseBed(.c, 
                             mingap = dist_cutoff, 
                             return_idx_only = T)
    .c[, c("start", "end"):= .(min(start), max(end)), idx]
    .c$idx <- NULL
    .c <- .c[, .SD[which.max(qValue)], .(seqnames, start, end)]
    fwrite(.c,
           filtered_peaks,
           sep= "\t",
           quote= F,
           na= ".",
           col.names= F)
  }
  print("DONE")
}, .(filtered_peaks, peaks_merge, dist_cutoff)]

# Merged_peaks ----
meta[!is.na(merged_file), {
  if(!file.exists(merged_file))
  {
    .c <- vl_importBed(unique(filtered_peaks))
    .c$idx <- vl_collapseBed(.c, 
                             mingap = dist_cutoff, 
                             return_idx_only = T)
    .c[, c("start", "end"):= .(min(start), max(end)), idx]
    .c$idx <- NULL
    .c <- .c[, .SD[which.max(qValue)], .(seqnames, start, end)]
    fwrite(.c,
           merged_file,
           sep= "\t",
           quote= F,
           na= ".",
           col.names= F)
  }
  print("DONE")
}, .(merged_file, dist_cutoff)]

# bw files ----
meta[!is.na(bw_file), {
  if(!file.exists(bw_file))
  {
    .g <- fread(gsub(ifelse(broad, "peaks.broadPeak$", "peaks.narrowPeak$"),
                     "treat_pileup.bdg" ,
                     peaks_rep),
                col.names = c("seqnames", "start", "end", "score"))
    .g[, start:= start+1]
    # Format GRanges
    .g <- GenomicRanges::GRanges(.g)
    GenomeInfoDb::seqlevels(.g, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
    GenomeInfoDb::seqlengths(.g) <- GenomeInfoDb::seqlengths(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
    # save
    rtracklayer::export.bw(.g,
                           con= bw_file)
  }
  print("done")
}, .(peaks_rep, bw_file, broad)]
meta[!is.na(bw_merge), {
  if(!file.exists(bw_merge))
  {
    vl_bw_merge(bw_file, 
                BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6, 
                bins_width = 50L, 
                output = bw_merge)
  }
  print("done")
}, bw_merge]

# PEAKS SAF files ----
meta[!is.na(peaks_saf), {
  if(any(!file.exists(peaks_saf)))
  {
    .c <- vl_importBed(merged_file)
    TSSs <- fread("db/saf/TSSs_0_0.saf", col.names = c("gene_id", "seqnames", "start", "end", "strand"))
    cl <- vl_closestBed(a = .c, b= TSSs)
    .c[cl, gene_id:= gene_id.b, on= c("seqnames", "start", "end")]
    .c <- .c[, .(GeneID= paste0("PEAK_", seqnames, ":", start, "-", end, ":", strand, "_", gene_id),
                 Chr= seqnames,
                 Start= start,
                 End= end,
                 Strand= strand)]
    fwrite(unique(.c), 
           peaks_saf, 
           sep= "\t")
  }
  print("done")
}, .(merged_file, peaks_saf, broad)]

# Compute counts ----
count_FUN <- function(saf, output, bam)
{
  saf <- fread(saf)
  saf <- saf[Chr %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")]
  saf <- as.data.frame(saf)
  .c <- featureCounts(bam, # count reads
                      annot.ext= saf,
                      isGTFAnnotationFile = F,
                      isPairedEnd = T,
                      nthreads = 8)
                      # allowMultiOverlap = F -> NOT USED, SHOULD BE CONSIDERED?
  saveRDS(.c, output)
}
meta[!is.na(peaks_counts), {
  if(!file.exists(peaks_counts))
    count_FUN(peaks_saf, peaks_counts, bam= bam)
  print("done")
}, .(peaks_saf, peaks_counts)]

# DESeq2 analysis
dds <- melt(meta, 
            id.vars = "ChIP", 
            measure.vars = patterns("counts$"), 
            variable.name = "feature", 
            value.name = "counts")
dds <- unique(na.omit(dds))
dds[, feature:= gsub("_counts$", "", feature)]
dds[, dds_file:= paste0("db/dds/cutnrun/", ChIP, "_", feature, ".dds")]
dds <- dds[, CJ(meta$cdition, "PH18", unique = T), (dds)]
dds <- dds[!V1 %in% c("PH18",  "WT") & ChIP %in% c("GAF", "PHO", "PC")]
dds[, FC_file:= paste0("db/FC_tables/cutnrun/", ChIP, "_", feature, "_", V1, "_vs_", V2, ".txt")]
dds[, {
  if(!file.exists(dds_file))
  {
    .c <- readRDS(counts)
    DF <- as.data.frame(.c[[1]])
    names(DF) <- gsub("[.]", "_", gsub(".bam$", "", names(DF)))
    DF <- DF[rowSums(DF)>100,]
    
    # Samples
    sampleTable <- as.data.frame(setNames(tstrsplit(names(DF), "_"), c("ChIP", "cdition", "rep")),
                                 row.names = names(DF))
    
    # DESeq2 dataset
    .dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                           colData= sampleTable,
                                           design= ~rep+cdition)
    # Libsize norm
    libsize <- .c[[4]][, -1]
    names(libsize) <- gsub("[.]", "_", gsub(".bam$", "", names(libsize)))
    libsize <- apply(libsize, 2, sum)
    DESeq2::sizeFactors(.dds) <- libsize/min(libsize)
    # Run DESeq2 and save
    .dds <- DESeq2::DESeq(.dds)
    saveRDS(.dds, dds_file)
  }else
    .dds <- readRDS(dds_file)
  
  # FC tables
  .SD[, {
    if(!file.exists(FC_file))
    {
      res <- as.data.frame(DESeq2::results(.dds,
                                           contrast= c("cdition", V1, V2)))
      res <- as.data.table(res, keep.rownames = "ID")
      fwrite(res,
             FC_file,
             col.names = T,
             sep= "\t",
             quote=F,
             na= NA)
    }
  }, .(V1, V2, FC_file)]
  print("DONE")
}, .(ChIP, counts, dds_file)]

# Add files to meta ----
add <- melt(dds, id.vars = c("ChIP", "feature", "V1"), measure.vars = patterns("file$"))
add[, name:= paste0(gsub("_file$", "", variable), "_", feature)]
add <- dcast(add, ChIP+V1~name, value.var = "value")
processed <- merge(meta, 
                   add, 
                   by.x= c("ChIP", "cdition"), 
                   by.y= c("ChIP", "V1"), 
                   all.x=T)

# SAVE ----
fwrite(processed,
       "Rdata/processed_metadata_CUTNRUN_review.txt",
       na= NA)

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(data.table)
require(parallel)
require(Rsubread)
require(seqinr)
require(rtracklayer)

#--------------------------------------------------------------#
# METDATA
#--------------------------------------------------------------#
meta <- readxl::read_xlsx("Rdata/metadata_cutnrun.xlsx")
meta <- as.data.table(meta)[!Comment %in% c("failed", "test")]
meta[is.na(Suffix), Suffix:= ""]

#--------------------------------------------------------------#
# Scer/Dm6 combined index
# Yeast genome -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna
#--------------------------------------------------------------#
if(length(list.files("/mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/", ".bt2$"))==0)
{
  dm6 <- read.fasta("/mnt/d/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa")
  S288C <- read.fasta("/mnt/d/_R_data/genomes/S288C/GCF_000146045.2_R64_genomic.fna")
  cmb <- c(dm6, S288C)
  write.fasta(sequences = cmb, 
              names = names(cmb),
              file.out = "/mnt/d/_R_data/genomes/dm6_S288C_combined/dm6_S288C_combined.fa")
  system("bowtie2-build /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C.fna /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C")
}

#--------------------------------------------------------------#
# Retrieve fastqs and trim reads
#--------------------------------------------------------------#
meta[, fq1:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = paste0(fq1, "$")), fq1]
meta[, fq2:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = paste0(fq2, "$")), fq2]

meta[, fq1_trim:= gsub(".fq.gz$", "_val_1.fq.gz", fq1)]
meta[, fq2_trim:= gsub(".fq.gz$", "_val_2.fq.gz", fq2)]

if(!all(file.exists(meta$fq1_trim, meta$fq2_trim)))
  mcmapply(FUN = function(outdir, fq1, fq2, fq1_trim, fq2_trim){
    if(!all(file.exists(c(fq1_trim, fq2_trim))))
      system(paste0("trim_galore --paired --gzip -o ", outdir, "/ ", fq1, " ", fq2))
    print("done")
  },
  outdir= dirname(meta$fq1),
  fq1= meta$fq1,
  fq2= meta$fq2,
  fq1_trim= meta$fq1_trim,
  fq2_trim= meta$fq2_trim,
  mc.preschedule = F,
  mc.cores = getDTthreads())

#--------------------------------------------------------------#
# Alignment
#--------------------------------------------------------------#
# Make chrom_sizes object
chrom_sizes <- rbind(data.table(fread("/mnt/d/_R_data/genomes/dm6/dm6.chrom.sizes.txt", 
                                      col.names = c("seqnames", "seqlengths")), 
                                type="ChIP"),
                     data.table(fread("/mnt/d/_R_data/genomes/S288C/S288C_contigs.txt"),
                                type="spike"))
setkeyv(chrom_sizes, "type")
meta[, bam:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/cutnrun/", ChIP, "_", cdition, "_", rep, Suffix, ".bam")]
meta[, {
  if(!file.exists(bam))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam)
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1_trim)
    cmd <- paste0(cmd, " -2 ", fq2_trim)
    cmd <- paste0(cmd, " -S ", sam_file)
    system(cmd)
    # sam to bam
    cmd <- paste0("/usr/bin/samtools view -@ 9 -b -q 30 ", sam_file, " -o ", bam)
    system(cmd)
    file.remove(sam_file)
  }
  print("DONE")
}, bam]

#--------------------------------------------------------------#
# Peak calling
#--------------------------------------------------------------#
meta[!ChIP %in% c("IgG", "input"), input:= ifelse(ChIP=="PH", "input", "IgG")]
meta[!is.na(input), input_bam:= meta[.BY, bam, on= c("ChIP==input", "rep", "cdition")], .(input, rep, cdition)]
MACS <- function(bam_ChIP,
                 bam_Input,
                 output,
                 broad,
                 SPMR)
{
  cmd <- paste0("/home/vloubiere/.local/bin/macs2 callpeak -g dm --keep-dup 1 -f BAMPE",
                ifelse(SPMR, " -B --SPMR", ""),
                " --outdir ", paste0(dirname(output), "/"),
                " -t ", paste(bam_ChIP, collapse= " "), 
                " -c ", paste(bam_Input, collapse= " "), 
                " -n ", basename(output))
  if(broad)
    cmd <- paste0(cmd, " --broad")
  return(cmd)
}
# Separated replicates
meta[!is.na(input), peaks_rep:= paste0('db/peaks/cutnrun/', paste0(ChIP, "_", cdition, "_", rep))]
meta[!is.na(input), broad:= ChIP %in% c("H3K27me3", "H2AK118Ub", "H3K4me1")]
meta[!is.na(input), {
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
# Merge
meta[!is.na(input), peaks_merge:= paste0('db/peaks/cutnrun/', paste0(ChIP, "_", cdition, "_merge"))]
meta[!is.na(input), {
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
# Return peaks file
meta[!is.na(input), peaks_rep:= paste0(peaks_rep, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]
meta[!is.na(input), peaks_merge:= paste0(peaks_merge, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]

#--------------------------------------------------------------#
# Confident peaks
#--------------------------------------------------------------#
meta[!is.na(input), filtered_peaks:= paste0("db/peaks/cutnrun/", ChIP, "_", cdition, 
                                            ifelse(broad, 
                                                   "_confident_peaks.broadPeak", 
                                                   "_confident_peaks.narrowPeak"))]
meta[ChIP=="H3K27me3", dist_cutoff:= 2500]
meta[ChIP=="H2AK118Ub", dist_cutoff:= 2500]
meta[ChIP=="H3K27Ac", dist_cutoff:= 250]
meta[ChIP=="H3K36me3", dist_cutoff:= 250]
meta[ChIP=="H3K4me1", dist_cutoff:= 250]
meta[ChIP=="PH", dist_cutoff:= 250]
meta[!is.na(input), {
  if(!file.exists(filtered_peaks))
  {
    .c <- vl_importBed(peaks_merge)
    .c <- .c[signalValue>2 & qValue>2]
    .c <- .c[vl_covBed(.c, unique(peaks_rep[rep=="rep1"]))>0]
    .c <- .c[vl_covBed(.c, unique(peaks_rep[rep=="rep2"]))>0]
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

#--------------------------------------------------------------#
# Merged_peaks
#--------------------------------------------------------------#
meta[!is.na(input), merged_file:= paste0("db/peaks/cutnrun_merged_peaks/", ChIP, "_merged_peaks", 
                                         ifelse(broad, ".broadPeak", ".narrowPeak"))]
meta[!is.na(input), {
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

#--------------------------------------------------------------#
# bw files
#--------------------------------------------------------------#
meta[!is.na(peaks_rep), bw_file:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_", rep, ".bw")]
meta[!is.na(peaks_rep), {
  if(!file.exists(bw_file))
  {
    .g <- fread(gsub(ifelse(broad, "peaks.broadPeak$", "peaks.narrowPeak$"),
                     "treat_pileup.bdg" ,
                     peaks_rep),
                col.names = c("seqnames", "start", "end", "score"))
    .g[, start:= start+1]
    # Format GRanges
    .g <- GRanges(.g)
    BS <- BSgenome::getBSgenome("dm6")
    GenomeInfoDb::seqlevels(.g, pruning.mode="coarse") <- GenomeInfoDb::seqlevels(BS)
    GenomeInfoDb::seqlengths(.g) <- GenomeInfoDb::seqlengths(BS)
    # save
    rtracklayer::export.bw(.g,
                           con= bw_file)
  }
  print("done")
}, .(peaks_rep, bw_file, broad)]
meta[!is.na(peaks_rep), bw_merge:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_merge.bw")]
meta[!is.na(peaks_rep), {
  if(!file.exists(bw_merge))
  {
    vl_bw_merge(bw_file, 
                "dm6", 
                bins_width = 25L, 
                output = bw_merge)
  }
  print("done")
}, .(peaks_rep, bw_file, broad)]

#--------------------------------------------------------------#
# Compute features SAF files
#--------------------------------------------------------------#
if(any(!file.exists(c("db/saf/promoters_750_250.saf",
                      "db/saf/geneBody_1000_end.saf",
                      "db/saf/TSSs_0_0.saf",
                      "db/saf/ATAC_peaks.saf"))))
{
  gtf <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
  seqlevelsStyle(gtf) <- "UCSC"
  gtf <- as.data.table(gtf)
  genes <- gtf[type=="gene" & gene_id %in% gtf[type %in% c("mRNA", "ncRNA"), gene_id]]
  # Promoters
  if(!file.exists("db/saf/promoters_750_250.saf"))
  {
    proms <- vl_resizeBed(genes, center = "start", upstream = 750, downstream = 250, genome = "dm6")
    proms <- proms[, .(GeneID= gene_id,
                       Chr= seqnames,
                       Start= start,
                       End= end,
                       Strand= strand)]
    fwrite(unique(proms), 
           "db/saf/promoters_750_250.saf", 
           sep= "\t")
  }
  # Gene body
  if(!file.exists("db/saf/geneBody_1000_end.saf"))
  {
    bodies <- vl_resizeBed(genes, center = "start", upstream = 1000, downstream = genes[, end-start+1], genome = "dm6")
    bodies <- bodies[, .(GeneID= gene_id,
                         Chr= seqnames,
                         Start= start,
                         End= end,
                         Strand= strand)]
    fwrite(unique(bodies), 
           "db/saf/geneBody_1000_end.saf", 
           sep= "\t")
  }
  # TTS
  if(!file.exists("db/saf/TTSs_2500_1000.saf"))
  {
    TTSs <- vl_resizeBed(genes, center = "end", upstream = 2500, downstream = 1000, genome= "dm6")
    TTSs <- TTSs[, .(GeneID= gene_id,
                     Chr= seqnames,
                     Start= start,
                     End= end,
                     Strand= strand)]
    fwrite(unique(TTSs), 
           "db/saf/TTSs_2500_1000.saf", 
           sep= "\t")
  }
  # TSS
  if(!file.exists("db/saf/TSSs_0_0.saf"))
  {
    TSSs <- vl_resizeBed(genes, center = "start", upstream = 0, downstream = 0, genome= "dm6")
    TSSs <- TSSs[, .(GeneID= gene_id,
                     Chr= seqnames,
                     Start= start,
                     End= end,
                     Strand= strand)]
    fwrite(unique(TSSs), 
           "db/saf/TSSs_0_0.saf", 
           sep= "\t")
  }
  # ATAC peaks
  if(!file.exists("db/saf/ATAC_peaks.saf"))
  {
    ATAC <- vl_importBed("db/peaks/ATAC/ATAC_confident_peaks.narrowPeak")
    TSSs <- fread("db/saf/TSSs_0_0.saf", col.names = c("gene_id", "seqnames", "start", "end", "strand"))
    cl <- vl_closestBed(a = ATAC, b= TSSs)
    ATAC[cl, gene_id:= gene_id.b, on= c("seqnames", "start", "end"), mult= "first"]
    ATAC <- ATAC[, .(GeneID= paste0("ATAC_", seqnames, ":", start, "-", end, ":", strand, "_", gene_id),
                     Chr= seqnames,
                     Start= start,
                     End= end,
                     Strand= strand)]
    fwrite(unique(ATAC), 
           "db/saf/ATAC_peaks.saf", 
           sep= "\t")
  }
}

#--------------------------------------------------------------#
# PEAKS SAF files
#--------------------------------------------------------------#
meta[!is.na(merged_file), peaks_saf:= paste0("db/saf/", ChIP, "_", "peaks.saf"), ChIP]
meta[!is.na(merged_file), {
  if(any(!file.exists(peaks_saf)))
  {
    .c <- vl_importBed(merged_file)
    TSSs <- fread("db/saf/TSSs_0_0.saf", col.names = c("gene_id", "seqnames", "start", "end", "strand"))
    cl <- vl_closestBed(a = .c, b= TSSs)
    .c[cl, gene_id:= gene_id.b, on= c("seqnames", "start", "end"), mult= "first"]
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

#--------------------------------------------------------------#
# Compute counts
#--------------------------------------------------------------#
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

meta[!is.na(input), peaks_counts:= paste0("db/counts/cutnrun/", ChIP, "_peaks_counts.txt")]
meta[!is.na(input), {
  if(!file.exists(peaks_counts))
    count_FUN(peaks_saf, peaks_counts, bam= bam)
  print("done")
}, .(peaks_saf, peaks_counts)]

meta[ChIP %in% c("PH", "H3K27Ac"), atac_counts:= paste0("db/counts/cutnrun/", ChIP, "_atac_counts.txt")]
meta[ChIP %in% c("PH", "H3K27Ac"), {
  if(!file.exists(atac_counts))
    count_FUN("db/saf/ATAC_peaks.saf", atac_counts, bam= bam)
  print("done")
}, atac_counts]

meta[ChIP %in% c("PH", "H3K27Ac"), prom_counts:= paste0("db/counts/cutnrun/", ChIP, "_prom_counts.txt")]
meta[ChIP %in% c("PH", "H3K27Ac"), {
  if(!file.exists(prom_counts))
    count_FUN("db/saf/promoters_750_250.saf", prom_counts, bam= bam)
  print("done")
}, prom_counts]

meta[ChIP %in% c("H3K27me3", "H2AK118Ub", "H3K4me1"), body_counts:= paste0("db/counts/cutnrun/", ChIP, "_body_counts.txt")]
meta[ChIP %in% c("H3K27me3", "H2AK118Ub", "H3K4me1"), {
  if(!file.exists(body_counts))
    count_FUN("db/saf/geneBody_1000_end.saf", body_counts, bam= bam)
  print("done")
}, body_counts]

meta[ChIP=="H3K36me3", TTS_counts:= paste0("db/counts/cutnrun/", ChIP, "_TTS_counts.txt")]
meta[ChIP=="H3K36me3", {
  if(!file.exists(TTS_counts))
    count_FUN("db/saf/TTSs_2500_1000.saf", TTS_counts, bam= bam)
  print("done")
}, TTS_counts]

#--------------------------------------------------------------#
# DESeq2 analysis
#--------------------------------------------------------------#
dds <- melt(meta, 
            id.vars = "ChIP", 
            measure.vars = patterns("counts$"), 
            variable.name = "feature", 
            value.name = "counts")
dds <- unique(na.omit(dds))
dds[, feature:= gsub("_counts$", "", feature)]
dds[, dds_file:= paste0("db/dds/cutnrun/", ChIP, "_", feature, ".dds")]
dds <- dds[, CJ(meta$cdition, "PH18", unique = T), (dds)]
dds <- dds[V1!="PH18"]
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
    sizeFactors(.dds) <- libsize/min(libsize)
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

#--------------------------------------------------------------#
# Add files to meta
#--------------------------------------------------------------#
add <- melt(dds, id.vars = c("ChIP", "feature", "V1"), measure.vars = patterns("file$"))
add[, name:= paste0(gsub("_file$", "", variable), "_", feature)]
add <- dcast(add, ChIP+V1~name, value.var = "value")
processed <- merge(meta, 
                   add, 
                   by.x= c("ChIP", "cdition"), 
                   by.y= c("ChIP", "V1"), 
                   all.x=T)

#--------------------------------------------------------------#
# SAVE
#--------------------------------------------------------------#
fwrite(processed,
       "Rdata/processed_metadata_CUTNRUN.txt",
       na= NA)


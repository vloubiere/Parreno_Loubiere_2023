setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(data.table)
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
# Alignment
#--------------------------------------------------------------#
# Make chrom_sizes object
chrom_sizes <- rbind(data.table(fread("/mnt/d/_R_data/genomes/dm6/dm6.chrom.sizes.txt", 
                                      col.names = c("seqnames", "seqlengths")), 
                                type="ChIP"),
                     data.table(fread("/mnt/d/_R_data/genomes/S288C/S288C_contigs.txt"),
                                type="spike"))
setkeyv(chrom_sizes, "type")
# Processing
meta[, fq1:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = fq1), fq1]
meta[, fq2:= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                        recursive = T, 
                        full.names = T, 
                        pattern = fq2), fq2]
meta[, bam:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/cutnrun/", ChIP, "_", cdition, "_", rep, Suffix, ".bam")]
meta[, {
  if(!file.exists(bam))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam)
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1)
    cmd <- paste0(cmd, " -2 ", fq2)
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
# Bed files
#--------------------------------------------------------------#
meta[, ChIP_bed:= paste0("db/bed/cutnrun/", ChIP, "_", cdition, "_", rep, Suffix, "_uniq.bed")]
meta[, spikein_bed:= paste0("db/bed/cutnrun/", ChIP, "_", cdition, "_", rep, Suffix, "_uniq_spikein.bed")]
meta[, {
  if(!file.exists(ChIP_bed))
  {
    # Import as bed
    bed <- fread(cmd= paste0("/usr/bin/samtools view -@ 9 -b -q 30 ", bam, " | bedtools bamtobed -i stdin -bedpe"))
    # Clean bed
    reads <- unique(bed[, .(seqnames= V1, start= V2, end= V6)])
    setorderv(reads, c("seqnames", "start", "end"))
    
    # Split and save
    vl_exportBed(reads[seqnames %in% chrom_sizes["ChIP", seqnames], .(seqnames, start, end, strand)], 
                 ChIP_bed)
    vl_exportBed(reads[seqnames %in% chrom_sizes["spike", seqnames], .(seqnames, start, end, strand)], 
                 spikein_bed)
  }
  print("DONE")
}, .(ChIP_bed, spikein_bed, bam)]

#--------------------------------------------------------------#
# bw files
#--------------------------------------------------------------#
meta[, bw:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_", rep, Suffix, ".bw")]
meta[, {
  if(!file.exists(bw))
  {
    .b <- GRanges(vl_importBed(ChIP_bed))
    cov <- GenomicRanges::coverage(.b)/length(.b)*1e6
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw)
  }
  print("DONE")
}, bw]
meta[, bw_merge:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, Suffix, "_merge.bw")]
meta[, {
  if(!file.exists(bw_merge))
  {
    .b <- GRanges(vl_importBed(ChIP_bed))
    cov <- GenomicRanges::coverage(.b)/length(.b)*1e6
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw_merge)
  }
  print("DONE")
}, bw_merge]

#--------------------------------------------------------------#
# Peak calling
#--------------------------------------------------------------#
meta[!ChIP %in% c("IgG", "input"), input:= ifelse(ChIP=="PH", "input", "IgG")]
meta[!is.na(input), input_bam:= meta[.BY, bam, on= c("ChIP==input", "rep", "cdition")], .(input, rep, cdition)]
MACS <- function(bam_ChIP,
                 bam_Input,
                 output,
                 broad)
{
  cmd <- paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 -g dm --keep-dup 1 -f BAMPE --outdir ", 
                paste0(dirname(output), "/"),
                " -t ", paste(bam_ChIP, collapse= " "), 
                " -c ", paste(bam_Input, collapse= " "), 
                " -n ", basename(output))
  if(broad)
    cmd <- paste0(cmd, " --broad")
  return(cmd)
}
# Separated replicates
meta[!is.na(input), peaks_rep:= paste0('db/peaks/cutnrun/', paste0(ChIP, "_", cdition, "_", rep))]
meta[!is.na(input), broad:= ChIP %in% c("H3K27me3", "H2AK118Ub", "H3K36me3", "H3K4me1")]
meta[!is.na(input), {
  if(!file.exists(paste0(peaks_rep, "_peaks.xls")))
  {
    cmd <- MACS(bam, 
                input_bam, 
                peaks_rep, 
                broad = broad)
    system(cmd)
    print("DONE")
  }
}, .(ChIP, peaks_rep, broad)]
meta[!is.na(input), peaks_rep:= paste0(peaks_rep, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]
# Merge
meta[!is.na(input), peaks_merge:= paste0('db/peaks/cutnrun/', paste0(ChIP, "_", cdition, "_merge"))]
meta[!is.na(input), {
  if(!file.exists(paste0(peaks_merge, "_peaks.xls")))
  {
    cmd <- MACS(bam, 
                input_bam, 
                peaks_merge, 
                broad = ChIP %in% c("H3K27me3", "H2AK118Ub", "H3K36me3", "H3K4me1"))
    system(cmd)
    print("DONE")
  }
}, .(ChIP, peaks_merge, broad)]
meta[!is.na(input), peaks_merge:= paste0(peaks_merge, ifelse(broad, "_peaks.broadPeak", "_peaks.narrowPeak"))]

#--------------------------------------------------------------#
# Confident peaks
#--------------------------------------------------------------#
meta[!is.na(input), filtered_peaks:= paste0("db/peaks/cutnrun/", ChIP, "_", cdition, "_confident_peaks.bed")]
meta[!is.na(input), {
  if(!file.exists(filtered_peaks))
  {
    .c <- vl_importBed(unique(peaks_merge))[,1:9]
    .c <- .c[vl_covBed(.c, unique(peaks_rep[rep=="rep1"]))>0]
    .c <- .c[vl_covBed(.c, unique(peaks_rep[rep=="rep2"]))>0]
    setnames(.c, paste0("V", seq(ncol(.c))))
    fwrite(.c,
           filtered_peaks,
           sep= "\t",
           quote= F,
           na= ".",
           col.names= F)
  }
  print("DONE")
}, filtered_peaks]

#--------------------------------------------------------------#
# Merged_peaks
#--------------------------------------------------------------#
meta[ChIP=="H3K27me3", c("enr_cutoff", "dist_cutoff"):= .(2, 2500)]
meta[ChIP=="H2AK118Ub", c("enr_cutoff", "dist_cutoff"):= .(2, 2500)]
meta[ChIP=="H3K27Ac", c("enr_cutoff", "dist_cutoff"):= .(3, 250)]
meta[ChIP=="H3K36me3", c("enr_cutoff", "dist_cutoff"):= .(3, 250)]
meta[ChIP=="H3K4me1", c("enr_cutoff", "dist_cutoff"):= .(3, 500)]
meta[ChIP=="PH", c("enr_cutoff", "dist_cutoff"):= .(3, 150)]
meta[!is.na(input), merged_file:= paste0("db/peaks/cutnrun_merged_peaks/", ChIP, "_merged_peaks.bed")]
meta[!is.na(input), {
  if(!file.exists(merged_file))
  {
    .c <- vl_importBed(filtered_peaks)
    .c <- vl_collapseBed(.c[V7>enr_cutoff & V9>2], 
                         mingap = dist_cutoff)
    vl_exportBed(.c, merged_file)
  }
  print("DONE")
}, .(merged_file, filtered_peaks, enr_cutoff, dist_cutoff)]

#--------------------------------------------------------------#
# Compute counts merged peaks
#--------------------------------------------------------------#
meta[!is.na(input), read_counts:= paste0("db/counts/cutnrun/", ChIP, "_counts.txt")]
meta[!is.na(input), {
  if(!file.exists(read_counts))
  {
    peaks <- vl_importBed(merged_file)
    files <- unique(ChIP_bed)
    names <- gsub("_uniq.bed", "", basename(files))
    peaks[, (names):= lapply(files, function(x) vl_covBed(peaks, x))]
    fwrite(peaks, 
           file = read_counts,
           quote= F,
           sep= "\t",
           col.names = T)
  }
  print("DONE")
}, .(merged_file, read_counts)]

#--------------------------------------------------------------#
# DESeq2 analysis
#--------------------------------------------------------------#
meta[!is.na(input), dds_file:= paste0("db/dds/cutnrun/", ChIP, ".dds"), ChIP]
meta[!is.na(input), {
  if(!file.exists(dds_file))
  {
    counts <- fread(read_counts)
    cols <- grep("_PH", names(counts))
    counts <- counts[rowSums(counts[, cols, with= F])>100]
    
    # Format
    DF <- data.frame(counts[, cols, with= F], 
                     row.names = counts[, paste0(seqnames, ":", start, "-", end)])
    sampleTable <- as.data.frame(setNames(tstrsplit(names(DF), "_"), c("ChIP", "cdition", "rep")),
                                 row.names = names(DF))
    
    # Run DESeq2 and save dds
    dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                          colData= sampleTable,
                                          design= ~rep+cdition)
    libsize <- sapply(colnames(dds), function(x) {
      cmd <- paste("wc -l", list.files("db/bed/cutnrun/", paste0(x, "_uniq.bed"), full.names = T))
      fread(cmd = cmd)$V1
    })
    sizeFactors(dds) <- libsize/min(libsize)
    dds <- DESeq2::DESeq(dds)
    saveRDS(dds, dds_file)
  }
  print("DONE")
}, .(ChIP, read_counts, dds_file)]

# FC tables
meta[!is.na(input) & cdition!="PH18", FC_file:= paste0("db/FC_tables/cutnrun/", ChIP, "_", cdition, "_vs_PH18.txt")]
meta[!is.na(FC_file),
     {
       if(!file.exists(FC_file))
       {
         dds <- readRDS(dds_file)
         res <- as.data.frame(DESeq2::results(dds, 
                                              contrast= c("cdition", cdition, "PH18")))
         res <- as.data.table(res, keep.rownames = T)
         res[, c("seqnames", "start", "end"):= tstrsplit(rn, ":|-")]
         fwrite(res[, .(seqnames, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)],
                FC_file,
                col.names = T, 
                sep= "\t",
                quote=F)
       }
       print("DONE")
     }, .(dds_file, FC_file, cdition)]

#--------------------------------------------------------------#
# SAVE
#--------------------------------------------------------------#
fwrite(meta, 
       "Rdata/processed_metadata_CUTNRUN.txt", 
       na= NA)


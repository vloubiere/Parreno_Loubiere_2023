setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

meta <- as.data.table(read_xlsx("Rdata/metadata_ATAC.xlsx"))
meta[, fq_file:= list.files("db/fastq/ATAC_JJ_2018/", fq_basename, full.names = T), fq_basename]
meta[, sam_file:= paste0("db/sam/ATAC/", gsub(".fastq.gz", ".sam", fq_basename)), fq_basename]
meta[, bed_file:= paste0("db/bed/ATAC/", gsub(".fastq.gz", ".bed", fq_basename)), fq_basename]
meta[, peaks_file:= paste0("db/narrowpeaks/ATAC/", gsub(".fastq.gz$", "", fq_basename), "_peaks.narrowPeak"), fq_basename]
meta[, bw_file:= paste0("db/bw/ATAC/", gsub(".fastq.gz", ".bw", fq_basename)), fq_basename]

meta[, {
  if(!file.exists(sam_file))
  {
    # bowtie 2
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6/Sequence/Bowtie2Index/genome --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -U ", fq_file)
    cmd <- paste0(cmd, " -S ", sam_file)
    system(cmd)
  }
  if(!file.exists(bed_file))
  {
    # Import as bed
    bed <- fread(cmd= paste0("/usr/local/bin/samtools view -@ 9 -b -q 30 ", sam_file, " | bedtools bamtobed -i stdin"))
    
    # Clean bed
    reads <- unique(bed[, .(seqnames= V1, start= ifelse(V6=="+", V2, V3-300), strand= V6)])
    reads[, end:= start+300]
    setkeyv(reads, c("seqnames", "start"))
    # Save
    vl_exportBed(reads, bed_file)
  }
  # Peak calling
  if(!file.exists(peaks_file))
  {
    cmd <- paste0("/home/vloubiere/.local/bin/macs3 callpeak -t ", bed_file)
    cmd <- paste0(cmd, " --keep-dup 1 -g dm --nomodel --outdir db/narrowpeaks/ATAC/ --name ")
    cmd <- paste0(cmd, gsub("_peaks.narrowPeak", "", basename(peaks_file)))
    system(cmd)
  }
  # Generate bw
  if(!file.exists(bw_file))
  {
    .b <- import(bed_file)
    cov <- GenomicRanges::coverage(.b)/length(.b)*1e6
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw_file)
  }
  print("DONE")
}, .(sam_file, bed_file, peaks_file, bw_file)]

if(!file.exists("db/narrowpeaks/ATAC/ATAC_merged_peaks.narrowPeak"))
{
  cmd <- paste0("/home/vloubiere/.local/bin/macs3 callpeak -t ", paste0(meta$bed_file, collapse= " "))
  cmd <- paste0(cmd, " -g dm --nomodel --outdir db/narrowpeaks/ATAC/ --name ATAC_merged")
  system(cmd)
  merged_peaks <- vl_importBed("db/narrowpeaks/ATAC/ATAC_merged_peaks.narrowPeak")
  rep_peaks <- rbindlist(lapply(meta$peaks_file, vl_importBed), idcol = T)
  merged_peaks$N_replicates <- rep_peaks[merged_peaks, length(unique(.id)), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
  vl_exportBed(merged_peaks[N_replicates>=12, 1:10][order(seqnames, start)],
               "db/narrowpeaks/ATAC/ATAC_filtered_merged_peaks.narrowPeak")
}

if(!file.exists("db/bw/ATAC/ATAC_merged.bw"))
{
  vl_bw_merge(x = meta$bw_file, 
              output = "db/bw/ATAC/ATAC_merged.bw", 
              BSgenome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6,
              bins = 50)
  }
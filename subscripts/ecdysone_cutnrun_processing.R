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
meta <- readxl::read_xlsx("Rdata/metadata_ecdysone_cutnrun.xlsx")
meta <- as.data.table(meta)
meta[, fastq:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/cutnrun_EcR/", fastq)]

#--------------------------------------------------------------#
# Download combined fastq files (contain read1+2)
#--------------------------------------------------------------#
meta[, {
  if(!file.exists(fastq))
    download.file(download_fq, 
                  fastq)
  print("done")
}, fastq]

#--------------------------------------------------------------#
# Split reads 1 and 2
#--------------------------------------------------------------#
meta[, c("fq1", "fq2"):= lapply(c("_1.fastq.gz", "_2.fastq.gz"), function(x) gsub(".fastq.gz$", x, fastq))]
meta[, {
  if(!file.exists(fq1))
  {
    .c <- fread(fastq, header= F, col.names = "fq")
    .c[, read:= rep(grepl(".1 ", grep("^@", fq, value = T), fixed = T), each=4)]
    .c[grepl("+SRR", fq, fixed = T), fq:= "+"]
    .c[grepl("^@", fq) & (read), fq:= gsub(".1 ", " ", fq, fixed = T)]
    .c[grepl("^@", fq) & !(read), fq:= gsub(".2 ", " ", fq, fixed = T)]
    tmp1 <- tempfile(fileext = ".fastq.gz")
    fwrite(.c[(read), .(fq)],
           tmp1,
           compress= "gzip",
           col.names = F)
    tmp2 <- tempfile(fileext = ".fastq.gz")
    fwrite(.c[!(read), .(fq)],
           tmp2,
           compress= "gzip",
           col.names = F)
    # Trim adapters
    cmd <- paste0("/usr/bin/fastp -i ", tmp1, 
                  " -I ", tmp2, 
                  " -o ", fq1, 
                  " -O ", fq2, " -g -x -p")
    system(cmd)
  }
}, fastq]

#--------------------------------------------------------------#
# Alignment
#--------------------------------------------------------------#
meta[, bam:= paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/cutnrun_EcR/EcR_", cdition, "_rep", rep, ".bam")]
meta[, {
  if(!file.exists(bam))
  {
    # bowtie 2
    sam_file <- gsub(".bam$", ".sam", bam)
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6/bowtie2_idx/BDGP6 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", paste0(fq1, collapse = ","))
    cmd <- paste0(cmd, " -2 ", paste0(fq2, collapse = ","))
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
meta[, ChIP_bed:= paste0("db/bed/cutnrun_EcR/EcR_", cdition, "_rep", rep, "_uniq.bed")]
meta[, {
  if(!file.exists(ChIP_bed))
  {
    # Import as bed
    bed <- fread(cmd= paste0("/usr/bin/samtools view -@ 9 -b -q 30 ", bam, " | bedtools bamtobed -i stdin -bedpe"))
    # Clean bed
    reads <- unique(bed[, .(seqnames= V1, start= V2, end= V6)])
    reads[, seqnames:= paste0("chr", seqnames)]
    setorderv(reads, c("seqnames", "start", "end"))
    
    # Split and save
    vl_exportBed(reads[, .(seqnames, start, end, strand)], 
                 ChIP_bed)
  }
  print("DONE")
}, .(ChIP_bed, bam)]

#--------------------------------------------------------------#
# bw files
#--------------------------------------------------------------#
meta[, bw:= paste0("db/bw/cutnrun_EcR/EcR_", cdition, "_rep", rep, ".bw")]
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
meta[, bw_merge:= paste0("db/bw/cutnrun_EcR/EcR_", cdition, "_merge.bw")]
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
# SAVE
#--------------------------------------------------------------#
fwrite(meta, 
       "Rdata/processed_metadata_ecdysone_cutnrun.txt", 
       na= NA)


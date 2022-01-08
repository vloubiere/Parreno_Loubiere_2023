setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(Rsubread)
require(seqinr)
require(rtracklayer)

#--------------------------------------------------------------#
# METDATA
#--------------------------------------------------------------#
meta <- data.table(file= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                                    ".fq.gz$", 
                                    full.names = T,
                                    recursive = T))
meta[grepl("^Ac", basename(file)), ChIP:= "H3K27Ac"]
meta[grepl("^me3", basename(file)), ChIP:= "H3K27me3"]
meta[grepl("18_|10_", basename(file)), cdition:= "PH18"]
meta[grepl("09_", basename(file)), cdition:= "PHD9"]
meta[grepl("11_", basename(file)), cdition:= "PHD11"]
meta[grepl("29_", basename(file)), cdition:= "PH29"]
meta[grepl("_1_", basename(file)), rep:= "rep1"]
meta[grepl("_2_", basename(file)), rep:= "rep2"]
meta[, sam:= paste0("db/sam/cutnrun/", ChIP, "_", cdition, "_", rep, ".sam")]
meta[, ChIP_bed:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/reps/", ChIP, "_", cdition, "_", rep, "_uniq.bed")]
meta[, spikein_bed:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/reps/", ChIP, "_", cdition, "_", rep, "_uniq_spikein.bed")]
meta[, ChIP_bed_merge:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/merge/", ChIP, "_", cdition, "_merge_uniq.bed")]
meta[, spikein_bed_merge:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/merge/", ChIP, "_", cdition, "_merge_uniq_spikein.bed")]
meta[, bw_reps:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bw/cutnrun_reps_vl/", ChIP, "_", cdition, "_", rep, ".bw"), ]
meta[, bw_merge:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bw/cutnrun_merge_vl/", ChIP, "_", cdition, "_merge.bw"), ]
fwrite(meta, 
       "Rdata/metadata_cutnrun_final.txt")

#--------------------------------------------------------------#
# Scer/Dm6 combined index
# Yeast genome -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.fna
#--------------------------------------------------------------#
Rsubread_index_prefix <- "/mnt/d/_R_data/genomes/dm6_S288C_combined/dm6_S288C_combined_index/dm6_S288C_combined_index"
if(F)
{
  dm6 <- read.fasta("/mnt/d/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa")
  S288C <- read.fasta("/mnt/d/_R_data/genomes/S288C/GCF_000146045.2_R64_genomic.fna")
  cmb <- c(dm6, S288C)
  write.fasta(sequences = cmb, 
              names = names(cmb),
              file.out = "/mnt/d/_R_data/genomes/dm6_S288C_combined/dm6_S288C_combined.fa")
  Rsubread::buildindex(basename = Rsubread_index_prefix, 
                       reference = "/mnt/d/_R_data/genomes/dm6_S288C_combined/dm6_S288C_combined.fa")
}

#--------------------------------------------------------------#
# Alignment
#--------------------------------------------------------------#
meta[, {
  if(!file.exists(sam))
  {
    Rsubread::align(index = Rsubread_index_prefix,
                    readfile1 = grep("_1.fq.gz$", file, value = T), 
                    readfile2 = grep("_2.fq.gz$", file, value = T), 
                    type= "dna", 
                    output_file = sam, 
                    maxMismatches= 6, 
                    maxFragLength = 1000,
                    unique= T, 
                    nTrim3 = 5,
                    nTrim5 = 5, 
                    nthreads= getDTthreads()-2, 
                    output_format = "SAM")
  }
}, sam]

#--------------------------------------------------------------#
# Filter aligned reads and genrate bed files
#--------------------------------------------------------------#
# Make chrom_sizes object
chrom_sizes <- rbind(cbind(fread("/mnt/d/_R_data/genomes/dm6/dm6.chrom.sizes.txt", 
                                 col.names = c("seqnames", "seqlengths")), 
                           data.table(type="ChIP")),
                     cbind(fread("/mnt/d/_R_data/genomes/S288C/S288C_contigs.txt"),
                           data.table(type="spike")))

# Clean reads and export bed split dm6/Scer
meta[, {
  if(!file.exists(ChIP_bed))
  {
    print("Import reads!")
    # Import
    reads <- data.table::fread(sam, 
                               fill= T, 
                               select = c("V1", "V2", "V3", "V4", "V5", "V10"), 
                               col.names = c("ID", "flag", "seqnames", "read_most_left_pos", "mapq", "read"))
    # Process
    reads[flag %in% c(83, 147), c("side", "pos"):= .("end", as.numeric(read_most_left_pos)+nchar(read))]
    reads[flag %in% c(99, 163), c("side", "pos"):= .("start", as.numeric(read_most_left_pos))]
    reads <- na.omit(reads)
    reads <- data.table::dcast(reads, ID+seqnames~side, value.var= list("pos", "mapq"))
    reads <- reads[mapq_start>=20 & mapq_end>=20, .(seqnames, start= pos_start, end= pos_end)]
    reads <- unique(reads)
    setorderv(reads, c("seqnames", "start"))
    # Split and save
    ChIP <- GenomicRanges::GRanges(reads[seqnames %in% chrom_sizes[type=="ChIP", seqnames], .(seqnames, start, end)])
    spike <- GenomicRanges::GRanges(reads[seqnames %in% chrom_sizes[type=="spike", seqnames], .(seqnames, start, end)])
    rtracklayer::export.bed(ChIP, ChIP_bed)
    rtracklayer::export.bed(spike, spikein_bed)
  }
  print("DONE")
}, .(ChIP_bed, spikein_bed, sam)]

#--------------------------------------------------------------#
# Merged bed files
#--------------------------------------------------------------#
meta[, {
  if(!file.exists(ChIP_bed_merge))
  {
    ChIP <- rbindlist(lapply(unique(ChIP_bed), fread))
    fwrite(ChIP, ChIP_bed_merge)
  }
  if(!file.exists(spikein_bed_merge))
  {
    spike <- rbindlist(lapply(unique(spikein_bed), fread))
    fwrite(spike, spikein_bed_merge)
  }
  print("DONE")
}, .(ChIP_bed_merge, spikein_bed_merge)]

#--------------------------------------------------------------#
# bw_files
#--------------------------------------------------------------#
meta[, {
  # Each rep separately
  .SD[, {
    if(!file.exists(bw_reps))
    {
      .b <- import(ChIP_bed)
      cov <- GenomicRanges::coverage(.b)/nrow(fread(spikein_bed))*1e4
      rtracklayer::export.bw(GRanges(cov), 
                             con= bw_reps)
    }
    print("DONE")
  }, .(ChIP_bed, spikein_bed, bw_reps)]
  # Merge
  if(!file.exists(bw_merge))
  {
    .b <- Reduce(c, lapply(unique(ChIP_bed), import))
    cov <- GenomicRanges::coverage(.b)/nrow(rbindlist(lapply(unique(spikein_bed), fread)))*1e4
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw_merge)
  }
  print("DONE")
}, bw_merge]


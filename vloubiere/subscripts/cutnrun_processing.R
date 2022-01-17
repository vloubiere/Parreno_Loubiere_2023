setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)
require(Rsubread)
require(seqinr)
require(rtracklayer)

#--------------------------------------------------------------#
# METDATA
#--------------------------------------------------------------#
meta <- data.table(fq1= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                                   "_1.fq.gz$", full.names = T,recursive = T),
                   fq2= list.files("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/Cut_n_run/", 
                                   "_1.fq.gz$", full.names = T,recursive = T))
meta[grepl("^Ac", basename(fq1)), ChIP:= "H3K27Ac"]
meta[grepl("^me3", basename(fq1)), ChIP:= "H3K27me3"]
meta[grepl("18_|10_", basename(fq1)), cdition:= "PH18"]
meta[grepl("09_", basename(fq1)), cdition:= "PHD9"]
meta[grepl("11_", basename(fq1)), cdition:= "PHD11"]
meta[grepl("29_", basename(fq1)), cdition:= "PH29"]
meta[grepl("_1_", basename(fq1)), rep:= "rep1"]
meta[grepl("_2_", basename(fq1)), rep:= "rep2"]
meta[, sam:= paste0("db/sam/cutnrun/", ChIP, "_", cdition, "_", rep, ".sam")]
meta[, ChIP_bed:= paste0("db/bed/cutnrun/reps/", ChIP, "_", cdition, "_", rep, "_uniq.bed")]
meta[, spikein_bed:= paste0("db/bed/cutnrun/reps/", ChIP, "_", cdition, "_", rep, "_uniq_spikein.bed")]
meta[, bw_reps:= paste0("db/bw/cutnrun_reps_vl/", ChIP, "_", cdition, "_", rep, ".bw")]
meta[, ChIP_bed_merge:= paste0("db/bed/cutnrun/merge/", ChIP, "_", cdition, "_merge_uniq.bed")]
meta[, spikein_bed_merge:= paste0("db/bed/cutnrun/merge/", ChIP, "_", cdition, "_merge_uniq_spikein.bed")]
meta[, bw_merge:= paste0("db/bw/cutnrun_merge_vl/", ChIP, "_", cdition, "_merge.bw"), ]
fwrite(meta, 
       "Rdata/processed_metadata_CUTNRUN.txt")

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
# Align and generate bed files
#--------------------------------------------------------------#
# Make chrom_sizes object
chrom_sizes <- rbind(data.table(fread("/mnt/d/_R_data/genomes/dm6/dm6.chrom.sizes.txt", 
                                      col.names = c("seqnames", "seqlengths")), 
                                type="ChIP"),
                     data.table(fread("/mnt/d/_R_data/genomes/S288C/S288C_contigs.txt"),
                                type="spike"))
setkeyv(chrom_sizes, "type")

# Processing
meta[, {
  if(!file.exists(sam))
  {
    # bowtie 2
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1)
    cmd <- paste0(cmd, " -2 ", fq2)
    cmd <- paste0(cmd, " -S ", sam)
    system(cmd)
  }
  if(!file.exists(ChIP_bed))
  {
    # Import as bed
    bed <- fread(cmd= paste0("/usr/local/bin/samtools view -@ 9 -b -q 30 ", sam, " | bedtools bamtobed -i stdin"))
    bed[, ID:= gsub("(.*)/.*", "\\1", V4)]
    # Clean bed
    reads <- merge(bed[V6=="+", .(seqnames= V1, start= V2, ID)], 
                   bed[V6=="-", .(seqnames= V1, end= V3, ID)], 
                   by= c("seqnames", "ID"))
    reads <- unique(reads[, .(seqnames, start, end)])
    setkeyv(reads, c("seqnames", "start"))
    # Split and save
    ChIP <- reads[seqnames %in% chrom_sizes["ChIP", seqnames], .(seqnames, start, end)]
    spike <- reads[seqnames %in% chrom_sizes["spike", seqnames], .(seqnames, start, end)]
    vl_exportBed(ChIP, ChIP_bed)
    vl_exportBed(spike, spikein_bed)
  }
  if(!file.exists(bw_reps))
  {
    .b <- import(ChIP_bed)
    .s <- import(spikein_bed)
    cov <- GenomicRanges::coverage(.b)/length(.s)*1e4
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw_reps)
  }
  print("DONE")
}, .(sam, ChIP_bed, spikein_bed, bw_reps)]

#--------------------------------------------------------------#
# Merged bed files
#--------------------------------------------------------------#
# Kick bad rep out of the merged
meta[!(ChIP=="H3K27me3" & rep=="rep2" & cdition=="PH18"), {
  if(!file.exists(ChIP_bed_merge))
  {
    ChIP <- vl_importBed(ChIP_bed)
    vl_exportBed(ChIP, ChIP_bed_merge)
  }
  if(!file.exists(spikein_bed_merge))
  {
    spike <- vl_importBed(spikein_bed)
    vl_exportBed(spike, spikein_bed_merge)
  }
  if(!file.exists(bw_merge))
  {
    ChIP <- import(ChIP_bed_merge)
    spike <- import(spikein_bed_merge)
    cov <- GenomicRanges::coverage(ChIP)/length(spike)*1e4
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw_merge)
  }
  print("DONE")
}, .(ChIP_bed_merge, spikein_bed_merge, bw_merge)]


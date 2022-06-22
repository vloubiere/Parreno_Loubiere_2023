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
meta <- as.data.table(meta)[Comment!="failed"]
meta[Suffix=="NA", Suffix:= ""]

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
    # browser()
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
    bed <- fread(cmd= paste0("/usr/bin/samtools view -@ 9 -b -q 30 ", bam, " | bedtools bamtobed -i stdin"))
    bed[, ID:= gsub("(.*)/.*", "\\1", V4)]
    # Clean bed
    reads <- bed[, .(seqnames= V1, 
                     start= min(V2), 
                     end= max(V3)), ID]
    reads <- unique(reads[, !"ID"])
    setkeyv(reads, c("seqnames", "start", "end"))
    
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
peak_calling <- merge(meta[ChIP %in% c("H3K27Ac", "H3K27me3", "H2AK118Ub") & Suffix=="", 
                           .(ChIP, cdition, rep, bam)],
                      meta[ChIP=="IgG" & Suffix=="", 
                           .(cdition, rep, bam)],
                      by= c("cdition", "rep"), 
                      allow.cartesian= T,
                      suffixes= c("_ChIP", "_Input"))
peak_calling[, output_prefix:= paste0(ChIP, "_", cdition, "_", rep)]
peak_calling[, cmd:= paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 -g dm --keep-dup 1 -f BAMPE --outdir ", 
                            normalizePath("db/peaks/K27_cutnrun/"),
                            " -t ", bam_ChIP, 
                            " -c ", bam_Input, 
                            " -n ", output_prefix)]
peak_calling[ChIP %in% c("H3K27me3", "H2AK118Ub"), cmd:= paste0(cmd, " --broad")]
#-------------------#
# Custom K118 
peak_calling[ChIP=="H2AK118Ub" & rep=="rep1" & cdition=="PH18", bam_ChIP:= "/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/cutnrun/H2AK118Ub_WT_rep1.bam"]
#-------------------#
peak_calling[, {
  check <- list.files("db/peaks/K27_cutnrun/", 
                      paste0(ChIP, "_", cdition, "_", rep, ".*peaks.xls"))
  if(length(check)<1)
    system(cmd)
  print("DONE")
}, cmd]
# Merge
peak_calling[, cmd:= paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 -g dm --keep-dup 1 -f BAMPE --outdir ", 
                            normalizePath("db/peaks/K27_cutnrun/"),
                            " -t ", paste(bam_ChIP, collapse= " "), 
                            " -c ", paste(bam_Input, collapse= " "),
                            " -n ", paste0(ChIP, "_", cdition, "_merge")), .(ChIP, cdition)]
peak_calling[ChIP %in% c("H3K27me3", "H2AK118Ub"), cmd:= paste0(cmd, " --broad")]
peak_calling[, {
  check <- list.files("db/peaks/K27_cutnrun/", 
                      paste0(ChIP, "_", cdition, "_merge_peaks.xls"))
  if(length(check)<1)
    system(cmd)
  print("DONE")
}, cmd]

# Confident peaks
peaks <- peak_calling[, .(file= list.files("db/peaks/K27_cutnrun/", 
                                           paste0(ChIP, "_", cdition, 
                                           fcase(ChIP=="H3K27me3", ".*.broadPeak$",
                                                 ChIP=="H2AK118Ub", ".*.broadPeak$",
                                                 ChIP=="H3K27Ac", ".*.narrowPeak$")),
                                           full.names = T)), .(ChIP, cdition)]
peaks[, rep:= fcase(grepl("rep1", file), "rep1",
                    grepl("rep2", file), "rep2",
                    grepl("merge", file), "merge")]
peaks[, filtered_peaks:= paste0("db/peaks/K27_cutnrun/", ChIP, "_", cdition, "_confident_peaks.bed")]
peaks[, {
  if(!file.exists(filtered_peaks))
  {
    .c <- fread(file[rep=="merge"])[, 1:9]
    .c$rep1 <- fread(file[rep=="rep1"])[.c, .N, .EACHI, on= c("V1", "V2<=V3", "V3>=V2")]$N
    .c$rep2 <- fread(file[rep=="rep2"])[.c, .N, .EACHI, on= c("V1", "V2<=V3", "V3>=V2")]$N
    fwrite(.c[rep1>0 & rep2>0, V1:V9],
           filtered_peaks,
           sep= "\t",
           quote= F)
  }
  print("DONE")
}, filtered_peaks]

# Merged_peaks
peaks[ChIP=="H3K27me3", c("enr_cutoff", "dist_cutoff"):= .(2, 2500)]
peaks[ChIP=="H2AK118Ub", c("enr_cutoff", "dist_cutoff"):= .(2, 2500)]
peaks[ChIP=="H3K27Ac", c("enr_cutoff", "dist_cutoff"):= .(3, 250)]
peaks[, merged_file:= paste0("db/peaks/K27_cutnrun_changes/", ChIP, "_merged_peaks.bed")]
peaks[, {
  .c <- vl_importBed(unique(filtered_peaks))
  .c <- vl_collapseBed(.c[V7>enr_cutoff & V9>2], 
                       mingap = dist_cutoff)
  vl_exportBed(.c, 
               merged_file)
}, .(merged_file, enr_cutoff, dist_cutoff)]

#--------------------------------------------------------------#
# Compute counts merged peaks
#--------------------------------------------------------------#
meta[ChIP %in% c("H3K27me3","H3K27Ac","H2AK118Ub") & Comment=="NA", read_counts:= paste0("db/counts/cutnrun/", ChIP, "_counts.txt")]

#-------------------#
# Custom K118 
meta[ChIP=="H2AK118Ub" & rep=="rep1" & cdition=="PH18", ChIP_bed:= "db/bed/cutnrun/H2AK118Ub_WT_rep1_uniq.bed"]
#-------------------#

meta[!is.na(read_counts), {
  if(!file.exists(read_counts))
  {
    peaks <- vl_importBed(paste0("db/peaks/K27_cutnrun_changes/", ChIP, "_merged_peaks.bed"))
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
}, .(ChIP, read_counts)]

#-------------------#
# Custom K118 
K118 <- fread("db/counts/cutnrun/H2AK118Ub_counts.txt")
setnames(K118, "H2AK118Ub_WT_rep1", "H2AK118Ub_PH18_rep1", skip_absent=TRUE)
fwrite(K118,
       "db/counts/cutnrun/H2AK118Ub_counts.txt",
       quote= F,
       sep= "\t",
       col.names = T)
#-------------------#

#--------------------------------------------------------------#
# DESeq2 analysis
#--------------------------------------------------------------#
meta[!is.na(read_counts), dds_file:= paste0("db/dds/cutnrun/", ChIP, ".dds"), ChIP]
meta[!is.na(read_counts) & ChIP=="H2AK118Ub", {
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
    #-------------------#
    # Custom K118
    if("H2AK118Ub_PH18_rep1" %in% names(libsize))
      libsize["H2AK118Ub_PH18_rep1"] <- fread(cmd = "wc -l db/bed/cutnrun/H2AK118Ub_WT_rep1_uniq.bed")$V1
    #-------------------#
    sizeFactors(dds) <- libsize/min(libsize)
    dds <- DESeq2::DESeq(dds)
    saveRDS(dds, dds_file)
  }
  dds <- readRDS(dds_file)
  for(temp in c("PH29", "PHD9", "PHD11"))
  {
    FC_file <- paste0("db/FC_tables/cutnrun/", ChIP, "_", temp, "_vs_PH18.txt")
    # if(!file.exists(FC_file))
    if(T)
    {
      res <- as.data.frame(DESeq2::results(dds, 
                                           contrast= c("cdition", temp, "PH18")))
      res <- as.data.table(res, keep.rownames = T)
      res[, c("seqnames", "start", "end"):= tstrsplit(rn, ":|-")]
      fwrite(res[, .(seqnames, start, end, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)],
             FC_file,
             col.names = T, 
             sep= "\t",
             quote=F)
    }
    print("DONE")
  }
}, .(ChIP, read_counts, dds_file)]
fwrite(meta, 
       "Rdata/processed_metadata_CUTNRUN.txt", 
       na= NA)


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
    cmd <- paste0("bowtie2 -p 10 -x /mnt/d/_R_data/genomes/dm6_S288C_combined_bowtie2/dm6_S288C --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700")
    cmd <- paste0(cmd, " -1 ", fq1)
    cmd <- paste0(cmd, " -2 ", fq2)
    cmd <- paste0(cmd, " -b ", bam)
    system(cmd)
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
meta[, bw:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_", rep, ".bw")]
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
meta[, bw_merge:= paste0("db/bw/cutnrun/", ChIP, "_", cdition, "_merge.bw")]
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
peak_calling <- merge(meta[ChIP %in% c("H3K27Ac", "H3K27me3"), 
                           .(ChIP,
                             cdition,
                             rep,
                             bam)],
                      meta[ChIP=="IgG" & Suffix=="", 
                           .(cdition, 
                             rep,
                             bam)],
                      by= c("cdition", "rep"), 
                      allow.cartesian= T,
                      suffixes= c("_ChIP", "_Input"))
peak_calling[, cmd:= paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 -g dm --keep-dup 1 -f BAMPE --outdir ", 
                            normalizePath("db/peaks/K27_cutnrun/"),
                            " -t ", bam_ChIP, 
                            " -c ", bam_Input, 
                            " -n ", paste0(ChIP, "_", cdition, "_", rep))]
peak_calling[ChIP=="H3K27me3", cmd:= paste0(cmd, " --broad")]
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
                            " -t ", paste(sam_ChIP, collapse= " "), 
                            " -c ", paste(sam_Input, collapse= " "),
                            " -n ", paste0(ChIP, "_", cdition, "_merge")), .(ChIP, cdition)]
peak_calling[ChIP=="H3K27me3", cmd:= paste0(cmd, " --broad")]
peak_calling[, {
  check <- list.files("db/peaks/K27_cutnrun/", 
                      paste0(ChIP, "_", cdition, "_merge_peaks.xls"))
  if(length(check)<1)
    system(cmd)
  print("DONE")
}, cmd]



peak_calling[, {
  #Peak calling parameters
  bs <- switch(ChIP, "H3K27Ac"= 100, "H3K27me3"= 1000)
  min_width <- switch(ChIP, "H3K27Ac"= 250, "H3K27me3"= 5000)
  min_dist <- switch(ChIP, "H3K27Ac"= 201, "H3K27me3"= 2501)
  min_qval <- switch(ChIP, "H3K27Ac"= 3, "H3K27me3"= 10)
  min_OR <- switch(ChIP, "H3K27Ac"= 1, "H3K27me3"= 2)

  # Peak calling
  if(!file.exists(peaks_narrowpeaks))
  {
    peaks <- vl_peakCalling(ChIP_bed_merge_ChIP,
                            ChIP_bed_merge_Input,
                            gaussian_blur = T,
                            BSgenome = BSgenome.Dmelanogaster.UCSC.dm6, 
                            bins_width = bs, 
                            bins_OR_cutoff = 1)
    # macs3 callpeak -t /mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed 
    # -c /mnt/d/_R_data/projects/epigenetic_cancer/db/bed/cutnrun/merge/IgG_PH18_merge_uniq.bed 
    # -g dm --nomodel --outdir /mnt/d/_R_data/projects/epigenetic_cancer/db/ --name test
    fwrite(peaks, 
           peaks_narrowpeaks, 
           col.names = F, 
           sep= "\t")
  }
  
  # Merged peaks
  if(!file.exists(peaks_merged))
  {
    peaks <- vl_importBed(peaks_narrowpeaks)
    # Merge (potentially over-segmented) peaks from peak calling
    .m <- vl_collapseBed(peaks, mingap = min_dist)[end-start>min_width]
    # Rec-omput enrichment of merged peaks
    merged <- vl_enrichBed(.m, 
                           ChIP_bed_merge, 
                           shuffled_bed)
    fwrite(merged[qValue>min_qval & signalValue>=min_OR], 
           peaks_merged,
           col.names = F, 
           sep= "\t")
  }
  print("done")
}, .(ChIP, cdition, ChIP_bed_merge_ChIP, ChIP_bed_merge_Input, peaks_narrowpeaks, peaks_merged)]

#--------------------------------------------------------------#
# Segment genome into homogenous regions and compute enrichment
#--------------------------------------------------------------#
meta[, dds_file:= paste0("db/dds/cutnrun/", ChIP, ".dds"), ChIP]
meta[, dds_file_merge:= paste0("db/dds/cutnrun/K27_segmentation_", ChIP, ".dds"), ChIP]

if(!file.exists("db/narrowpeaks/cutnrun/K27_segmentation_merged.bed"))
{
  if(!file.exists("db/narrowpeaks/cutnrun/K27Ac_mergedBed_peaks.narrowPeak"))
  {
    K27Ac_merged <- vl_peakCalling(ChIP= c("db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed",
                                           "db/bed/cutnrun/merge/H3K27Ac_PHD11_merge_uniq.bed",
                                           "db/bed/cutnrun/merge/H3K27Ac_PHD9_merge_uniq.bed",
                                           "db/bed/cutnrun/merge/H3K27Ac_PH29_merge_uniq.bed"),
                                   Input = c("db/bed/cutnrun/merge/IgG_PH18_merge_uniq.bed",
                                             "db/bed/cutnrun/merge/IgG_PHD11_merge_uniq.bed",
                                             "db/bed/cutnrun/merge/IgG_PHD9_merge_uniq.bed",
                                             "db/bed/cutnrun/merge/IgG_PH29_merge_uniq.bed"),
                                   bins_width = 100,
                                   BSgenome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
    K27Ac_merged <- vl_collapseBed(K27Ac_merged, mingap = 201)[end-start>101]
    vl_exportBed(K27Ac_merged, "db/narrowpeaks/cutnrun/K27Ac_mergedBed_peaks.narrowPeak")
  }else
    K27Ac_merged <- vl_importBed("db/narrowpeaks/cutnrun/K27Ac_mergedBed_peaks.narrowPeak")
  
  if(!file.exists("db/narrowpeaks/cutnrun/K27me3_mergedBed_peaks.narrowPeak"))
  {
    K27me3_merged <- vl_peakCalling(ChIP= c("db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed",
                                           "db/bed/cutnrun/merge/H3K27me3_PHD11_merge_uniq.bed",
                                           "db/bed/cutnrun/merge/H3K27me3_PHD9_merge_uniq.bed",
                                           "db/bed/cutnrun/merge/H3K27me3_PH29_merge_uniq.bed"),
                                   Input = c("db/bed/cutnrun/merge/IgG_PH18_merge_uniq.bed",
                                             "db/bed/cutnrun/merge/IgG_PHD11_merge_uniq.bed",
                                             "db/bed/cutnrun/merge/IgG_PHD9_merge_uniq.bed",
                                             "db/bed/cutnrun/merge/IgG_PH29_merge_uniq.bed"),
                                   bins_width = 1000,
                                   BSgenome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
    K27me3_merged <- vl_collapseBed(K27me3_merged, mingap = 2501)[end-start>5000]
    vl_exportBed(K27me3_merged, "db/narrowpeaks/cutnrun/K27me3_mergedBed_peaks.narrowPeak")
  }else
    K27me3_merged <- vl_importBed("db/narrowpeaks/cutnrun/K27me3_mergedBed_peaks.narrowPeak")
  
  # Collapse
  merged_peaks <- rbind(K27Ac_merged, K27me3_merged)
  merged_peaks <- vl_collapseBed(merged_peaks, return_idx_only = T)
  
  # Segment
  setkeyv(merged_peaks, c("seqnames", "start"))
  K27 <- merged_peaks[, {
    ranges <- sort(unique(c(start, end[-.N]+1, end[.N])))
    .(start= ranges[-length(ranges)], end= ranges[-1]-1)
  }, .(seqnames, idx)]
  K27 <- unique(K27)[end-start>300]
  
  # Compute WT enrichment
  K27$K27me3_log2_enr_PH18 <- vl_enrichBed(K27,
                                           ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed",
                                           Input_bed = "db/bed/cutnrun/merge/IgG_PH18_merge_uniq.bed")$signalValue
  K27$K27me3_log2_enr_PHD11 <- vl_enrichBed(K27,
                                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PHD11_merge_uniq.bed",
                                            Input_bed = "db/bed/cutnrun/merge/IgG_PHD11_merge_uniq.bed")$signalValue
  K27$K27me3_log2_enr_PHD9 <- vl_enrichBed(K27,
                                           ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PHD9_merge_uniq.bed",
                                           Input_bed = "db/bed/cutnrun/merge/IgG_PHD9_merge_uniq.bed")$signalValue
  K27$K27me3_log2_enr_PH29 <- vl_enrichBed(K27,
                                           ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PH29_merge_uniq.bed",
                                           Input_bed = "db/bed/cutnrun/merge/IgG_PH29_merge_uniq.bed")$signalValue
  K27$K27Ac_log2_enr_PH18 <- vl_enrichBed(K27,
                                          ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed",
                                          Input_bed = "db/bed/cutnrun/merge/IgG_PH18_merge_uniq.bed")$signalValue
  K27$K27Ac_log2_enr_PHD11 <- vl_enrichBed(K27,
                                           ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PHD11_merge_uniq.bed",
                                           Input_bed = "db/bed/cutnrun/merge/IgG_PHD11_merge_uniq.bed")$signalValue
  K27$K27Ac_log2_enr_PHD9 <- vl_enrichBed(K27,
                                          ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PHD9_merge_uniq.bed",
                                          Input_bed = "db/bed/cutnrun/merge/IgG_PHD9_merge_uniq.bed")$signalValue
  K27$K27Ac_log2_enr_PH29 <- vl_enrichBed(K27,
                                          ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PH29_merge_uniq.bed",
                                          Input_bed = "db/bed/cutnrun/merge/IgG_PH29_merge_uniq.bed")$signalValue
  
  # SAVE
  fwrite(K27[, !"idx"],
         col.names = T, 
         sep= "\t",
         "db/narrowpeaks/cutnrun/K27_segmentation_mergedBed.bed")
}

#--------------------------------------------------------------#
# DESeq2 analysis separate ChIP
#--------------------------------------------------------------#
# Make dds
FC <- meta[ChIP %in% c("H3K27me3","H3K27Ac"), {
  if(!file.exists(dds_file))
  {
    # Collapse peaks/ChIP
    merged_peaks <- vl_importBed(unique(peaks_merged))
    merged_peaks <- vl_collapseBed(merged_peaks)
    
    # Compute overlapping reads for all replicates and make DF
    .c <- lapply(unique(ChIP_bed), function(x) vl_covBed(merged_peaks, x))
    names(.c) <- paste0(ChIP, "_", cdition, "_", rep)
    DF <- data.frame(do.call(cbind, .c))
    rownames(DF) <- merged_peaks[, paste0(seqnames, ":", start, "-", end)]
    
    # Assemble DESeq2 SampleTable
    sampleTable <- data.table(colnames(DF), 
                              as.data.table(tstrsplit(colnames(DF), "_")))
    names(sampleTable) <- c("sample", "ChIP", "cdition", "rep")
    sampleTable <- data.frame(sampleTable[, .(cdition, rep)], 
                              row.names = sampleTable$sample)
    
    # Run DESeq2 and save dds
    dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                          colData= sampleTable,
                                          design= ~rep+cdition)
    dds <- DESeq2::DESeq(dds)
    libsize <- data.table(file= unique(ChIP_bed))[, fread(cmd= paste0("wc -l ", file)), file]$V1
    sizeFactors(dds) <- libsize/min(libsize)
    saveRDS(dds, dds_file)
    print("DONE")
  }else
    dds <- readRDS(dds_file)
  CJ(as.character(dds$cdition),
     as.character(dds$cdition), 
     unique= T)[V1!=V2]
}, .(ChIP, dds_file)]

FC <- FC[V2=="PH18"]
FC[, FC_file:= paste0("db/FC_tables/cutnrun/", ChIP, "_", V1, "_vs_", V2, ".txt"), .(ChIP, V1, V2)]
FC[, {
  dds <- readRDS(dds_file)
  .SD[, {
    if(!file.exists(FC_file))
    {
      # Compute FC tables
      .c <- DESeq2::lfcShrink(dds,
                              type= "ashr",
                              contrast= c("cdition", V1, V2))
      .c <- as.data.table(as.data.frame(.c), 
                          keep.rownames = "coor")
      fwrite(.c,
             FC_file, 
             col.names = T, 
             sep= "\t")
    }
    print("DONE")
  }, .(V1, V2, FC_file)]
  print("DONE")
}, dds_file]
meta[FC, FC_file_ChIP_peaks:= .(i.FC_file), on= c("ChIP", "cdition==V1")]

#--------------------------------------------------------------#
# DESeq2 analysis merged ChIP
#--------------------------------------------------------------#
FC_merge <- meta[ChIP %in% c("H3K27me3","H3K27Ac"), {
  if(!file.exists(dds_file_merge))
  {
    K27 <- vl_importBed("db/narrowpeaks/cutnrun/K27_segmentation_mergedBed.bed")[, 1:3]
      
    # Compute overlapping reads for all ChIP replicates
    .c <- lapply(unique(ChIP_bed), function(x) vl_covBed(K27, x))
    names(.c) <- paste0(ChIP, "_", cdition, "_", rep)
    
    # Assemble DESeq2 DF
    DF <- data.frame(do.call(cbind, .c))
    rownames(DF) <- K27[, paste0(seqnames, ":", start, "-", end)]
    
    # Assemble DESeq2 SampleTable
    sampleTable <- data.table(colnames(DF), 
                              as.data.table(tstrsplit(colnames(DF), "_")))
    names(sampleTable) <- c("sample", "ChIP", "cdition", "rep")
    sampleTable <- data.frame(sampleTable[, .(cdition, rep)], 
                              row.names = sampleTable$sample)
    
    # Run DESeq2 and save dds
    dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                          colData= sampleTable,
                                          design= ~rep+cdition)
    dds <- DESeq2::DESeq(dds)
    libsize <- data.table(file= unique(ChIP_bed))[, fread(cmd= paste0("wc -l ", file)), file]$V1
    sizeFactors(dds) <- libsize/min(libsize)
    saveRDS(dds, dds_file_merge)
    print("DONE")
  }else
    dds <- readRDS(dds_file_merge)
  CJ(as.character(dds$cdition),
     as.character(dds$cdition), 
     unique= T)[V1!=V2]
}, .(ChIP, dds_file_merge)]

FC_merge <- FC_merge[V2=="PH18"]
FC_merge[, FC_merge_file:= paste0("db/FC_tables/cutnrun/K27_segmentation_", ChIP, "_", V1, "_vs_", V2, ".txt"), .(ChIP, V1, V2)]
FC_merge[, {
  dds <- readRDS(dds_file_merge)
  .SD[, {
    if(!file.exists(FC_merge_file))
    {
      # Compute FC_merge tables
      .c <- DESeq2::lfcShrink(dds,
                              type= "ashr",
                              contrast= c("cdition", V1, V2))
      .c <- as.data.table(as.data.frame(.c), 
                          keep.rownames = "coor")
      fwrite(.c,
             FC_merge_file, 
             col.names = T, 
             sep= "\t")
    }
    print("DONE")
  }, .(V1, V2, FC_merge_file)]
  print("DONE")
}, dds_file_merge]


meta[FC_merge, FC_file_K27_segmentation:= .(i.FC_merge_file), on= c("ChIP", "cdition==V1")]
fwrite(meta, 
       "Rdata/processed_metadata_CUTNRUN.txt", 
       na = NA)

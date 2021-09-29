require(data.table)
require(vlfunctions)
require(BSgenome.Dmelanogaster.UCSC.dm6)

dat <- fread("Rdata/raw_metadata_final.txt")[cdition=="ATAC_ED" | grepl("H3K9me3", cdition)]

#--------------------------------------------#
# K9me3 and ATAC-Seq processing
#--------------------------------------------#
dir.create("D:/_R_data/projects/public_data/dm6/bw/ATAC_rep/", showWarnings = F)
dir.create("D:/_R_data/projects/public_data/dm6/bw/K9me3_rep/", showWarnings = F)

dat[, {
 bw_folder <- ifelse(grepl("ATAC", fq_file),
                     "D:/_R_data/projects/public_data/dm6/bw/ATAC_rep/",
                     "D:/_R_data/projects/public_data/dm6/bw/K9me3_rep/")
 .bn <- gsub(".fq$|.fastq$|.fq.gz$|.fastq.gz$", "", basename(fq_file))
 bw_ouptut <- paste0(bw_folder, .bn, '.bw')
 if(!file.exists(bw_ouptut))
 {
  vl_ChIP_pipeline(fq1 = fq_file,
                   chrom_sizes = fread("D:/_R_data/genomes/dm6/dm6.chrom.sizes_CHIP.txt",
                                       col.names = c("seqnames", "seqlengths")),
                   basename_prefix = "",
                   Rsubread_index_prefix = "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index",
                   bam_folder = "D:/_R_data/projects/public_data/dm6/bam/",
                   bw_folder = bw_folder,
                   extend = 300,
                   use_samtools = F)
  print(paste(bw_ouptut, "DONE!"))
 }else
  print(paste(bw_ouptut, "ALREADY EXISTS!"))
}, .(fq_file, prefix_output)]

#--------------------------------------------#
# Merged bw
#--------------------------------------------#
dir.create("D:/_R_data/projects/public_data/dm6/bw/merged_tracks", 
           showWarnings = F)
dat[, {
 output <- paste0("D:/_R_data/projects/public_data/dm6/bw/merged_tracks/", cdition, "_merged.bw")
 if(!file.exists(output))
 {
  files <- list.files("D:/_R_data/projects/public_data/dm6/bam/", 
                      paste0("^", prefix_output, ".*.bed$", collapse= "|"),
                      full.names = T)
  .c <- lapply(files, function(x) fread(x, 
                                        select = 1:3))
  .c <- rbindlist(.c)
  colnames(.c) <- c("seqnames", "start", "end")
  reads <- GRanges(.c)
  total_reads <- length(reads)
  cov <- GenomicRanges::coverage(reads)/total_reads * 1e+06
  
  rtracklayer::export.bw(GRanges(cov),
                         con = output)
  print(paste(output, "-->DONE"))
 }else
  print(paste(output, "-->ALREADY EXISTS"))
}, cdition]





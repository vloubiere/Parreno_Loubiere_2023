setwd("/_R_data/projects/epigenetic_cancer/")
require(Rsubread)
require(DESeq2)
require(data.table)

#----------------------------------------------------------#
# Build dm6 index
#----------------------------------------------------------#
# setwd("/_R_data/genomes/dm6/subreadr_index/")
# ref <- "/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
# buildindex(basename= "subreadr_dm6_index", reference= ref)
# setwd("/_R_data/projects/epigenetic_cancer/")

#----------------------------------------------------------#
# Metadata
#----------------------------------------------------------#
dat <- data.table(file= list.files("db/fastq", ".fq.gz", full.names = T, recursive = T))
dat[, cdition:= .(strsplit(sub("(_)(?=[^_]+$)", " ", .BY, perl=T), " ")[[1]][1]), gsub("rep|_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
dat[, replicate:= .(paste0("rep", strsplit(sub("(_)(?=[^_]+$)", " ", .BY, perl=T), " ")[[1]][2])), gsub("rep|_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
dat[, bam:= paste0("db/bam/", cdition, "_", replicate, ".bam"), .(cdition, replicate)]
dat[, counts_file:= paste0("db/counts/", gsub(".bam$", "_counts.rds", basename(bam)))]
dat[, group:= tstrsplit(file, "/", keep= 3)]
dat[, dds_file:= paste0("db/dds/", group, "_dds.rds")]

#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
dat[, {
  if(!file.exists(bam)){
    print("START...")
    if(.N==2){
      stats <- capture.output(subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", readfile1= file[1], readfile2= file[2], 
                                      maxMismatches = 6, nthreads = 10, unique = T, output_file= bam))
    }else if(.N==1){
      stats <- capture.output(subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", readfile1= file[1], 
                                      nthreads = 10, unique = T, output_file= bam))
    }
    writeLines(stats, con = gsub(".bam$", "_stats.txt", bam))
  }
  print(paste(bam, "DONE!"))
}, .(cdition, replicate, bam)]

#----------------------------------------------------------#
# Compute counts
#----------------------------------------------------------#
dat[, {
  if(!file.exists(counts_file)){
    if(.N==2){
      counts <- featureCounts(bam, annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", isGTFAnnotationFile = T, isPairedEnd = T, nthreads = 8)
    }else if(.N==1){
      counts <- featureCounts(bam, annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", isGTFAnnotationFile = T, isPairedEnd = F, nthreads = 8)
    }
    saveRDS(counts, counts_file)
  }
  print(paste(counts_file, "DONE!"))
}, .(cdition, replicate, bam, counts_file)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dat[, {
  sampleTable <- data.frame(unique(.SD[, .(cdition, replicate, counts_file= basename(counts_file))]), row.names = "counts_file")
  if(!file.exists(dds_file)){
    # DESEq2
    DF <- .SD[, {
      .c <- data.table(readRDS(counts_file)$counts, keep.rownames = T)
      colnames(.c)[2] <- "counts"
      .c
    }, counts_file]
    DF <- data.frame(dcast(DF, rn~basename(counts_file), value.var = "counts"), row.names = "rn")
    DF <- DF[rowSums(DF)>10, ]
    dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~replicate+cdition)
    dds <- DESeq(dds)
    saveRDS(dds, dds_file)
  }else{
    dds <- readRDS(dds_file)
  }
  print(paste(dds_file, "DONE!"))
  # Differential expression
  diff <- CJ(sampleTable$cdition, sampleTable$cdition, unique = T)
  diff <- diff[V1!=V2]
  diff[, FC_file:= paste0("db/FC_tables/", group, "_", V1, "_vs_", V2, "_FC.txt"), .(V1, V2)]
  diff[, {
    if(!file.exists(FC_file)){
      .c <- as.data.frame(lfcShrink(dds, type= "ashr", contrast= c("cdition", V1, V2)))
      # .c <- as.data.frame(results(dds, contrast= c("cdition", V1, V2)))
      fwrite(.c, FC_file, col.names = T, row.names = T, sep ="\t", quote= F)
    }
    print(paste(FC_file, "DONE!"))
  }, .(FC_file, V1, V2)]
}, .(group, dds_file)]








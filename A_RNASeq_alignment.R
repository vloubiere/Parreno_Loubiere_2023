setwd("/home/vloubiere/genomes/dm6/subreadr_index/")
require(Rsubread)

#----------------------------------------------------------#
# Build dm6 index
#----------------------------------------------------------#
# ref <- "/home/vloubiere/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
# buildindex(basename="subreadr_dm6_index", reference= ref)

#----------------------------------------------------------#
# Metadata
#----------------------------------------------------------#
dat <- data.table(file= list.files("/home/vloubiere/projects/epigenetic_cancer/db/fastq", ".fq.gz", full.names = T, recursive = T))
dat[, cdition:= .(strsplit(sub("(_)(?=[^_]+$)", " ", .BY, perl=T), " ")[[1]][1]), gsub("rep|_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
dat[, replicate:= .(paste0("rep", strsplit(sub("(_)(?=[^_]+$)", " ", .BY, perl=T), " ")[[1]][2])), gsub("rep|_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
dat[, bam:= paste0("/home/vloubiere/projects/epigenetic_cancer/db/bam/", cdition, "_", replicate, ".bam"), .(cdition, replicate)]
dat[, counts_file:= paste0("/home/vloubiere/projects/epigenetic_cancer/db/counts/", gsub(".bam$", "_counts.rds", basename(bam)))]
dat[, group:= tstrsplit(dat$file, "/", keep= 8)]

#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
dat[, {
  if(!file.exists(bam)){
    if(.N==2){
      stats <- capture.output(subjunc(index= "subreadr_dm6_index", readfile1= file[1], readfile2= file[2], 
                                      maxMismatches = 6, nthreads = 10, unique = T, output_file= bam))
    }else if(.N==1){
      stats <- capture.output(subjunc(index= "subreadr_dm6_index", readfile1= file[1], 
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
      counts <- featureCounts(bam, annot.ext= "/home/vloubiere/genomes/dm6/dmel-all-r6.36.gtf", isGTFAnnotationFile = T, isPairedEnd = T)
    }else if(.N==1){
      counts <- featureCounts(bam, annot.ext= "/home/vloubiere/genomes/dm6/dmel-all-r6.36.gtf", isGTFAnnotationFile = T, isPairedEnd = F)
    }
    saveRDS(counts, counts_file)
  }
  print(paste(counts_file, "DONE!"))
}, .(cdition, replicate, bam, counts_file)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dat[, {
  if(!file.exists(output)){
    sampleTable <- data.frame(unique(.SD[, .(cdition, replicate, counts_file)]), row.names = "counts_file")
    mat <- readRDS("Rdata/processed_peSTARRSeq_data/filtered_counts_prior_DESeq2.rds")
    sampleTable <- grep("rep", colnames(mat), value = T)
    sampleTable <- data.frame(condition= sapply(sampleTable, function(x) strsplit(x, "_")[[1]][1]),
                              replicate= sapply(sampleTable, function(x) strsplit(x, "_")[[1]][2]),
                              row.names= sampleTable)
    DF <- data.frame(mat[, DSCP_rep1:input_rep5], row.names = mat$rn)+1
    # DESeq2
    dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~replicate+condition)
    # SizeFactors
    sizeFactors(dds) <- estimateSizeFactorsForMatrix(as.matrix(DF[grep("control.*vs.*control", rownames(DF)),]))
    # Result
    res <- DESeq(dds)
    saveRDS(res, "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/dds_result_object.rds")
    
    # Differential expression
    diff <- as.data.table(as.data.frame(results(res, contrast= c("condition", "DSCP", "input"))), keep.rownames= T)
    diff[, c("enh_L", "enh_R"):= tstrsplit(rn, "_vs_")]
    diff <- diff[, .(enh_L, enh_R, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)]
    
    boxplot(diff[grepl("control", enh_L) & grepl("control", enh_R), log2FoldChange], notch= T)
    abline(h= 0, lty= 2)
    saveRDS(diff, "/groups/stark/vloubiere/projects/pe_STARRSeq/Rdata/processed_peSTARRSeq_data/DESeq2_FC_table.rds")
    
  }
  print(paste(output, "DONE!"))
}, .(group, output= paste0("/home/vloubiere/projects/epigenetic_cancer/db/dds/", group, "_dds.rds"))]








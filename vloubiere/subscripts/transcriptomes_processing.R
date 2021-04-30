#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
dat <- data.table(file= list.files("db/fastq", ".fq.gz", full.names = T, recursive = T))
dat[, project:= tstrsplit(file, "/", keep= 3)]
dat[, bam:= paste0("db/bam/", project, "/", .BY[1], ".bam"), gsub("_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
dat[, paired:= ifelse(.N==2, T, F), bam]
# Create directories
dat[, {
  .c <- paste0("db/bam/", project)
  if(!dir.exists(.c))
    dir.create(.c)
  }, project]
# align
dat[, {
  if(!file.exists(bam)) 
  {
    print("START...")
    subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", 
            readfile1= file[1], 
            readfile2= if(paired) file[2] else NULL, 
            maxMismatches = 6, 
            nthreads = 10, 
            unique = T, 
            output_file= bam)
  }
  print(paste(bam, "DONE!"))
}, .(bam, paired)]

#----------------------------------------------------------#
# Counts
#----------------------------------------------------------#
dat[, {
  .c <- paste0("db/counts/", project)
  if(!dir.exists(.c))
    dir.create(.c)
  }, project]
dat[, counts_file:= paste0("db/counts/", project, "/", gsub(".bam", "_counts.rds", basename(bam)))]
dat[, {
  if(!file.exists(counts_file))
  {
    counts <- featureCounts(bam, 
                            annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", 
                            isGTFAnnotationFile = T, 
                            isPairedEnd = ifelse(paired, T, F), 
                            nthreads = 8)
    saveRDS(counts, counts_file)
  }
  print(paste(counts_file, "DONE!"))
}, .(bam, counts_file, paired)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dat[, dds_file:= paste0("db/dds/", project, "_dds.rds")]
dat <- dat[!grepl("_T", file)] # Remove transplant experiments (NO REP)
dat[, rep:= gsub(".*_(.*)_counts.rds$", "\\1", counts_file)]
dat[, rep:= gsub("rep", "", rep)]
dat[, cdition:= gsub("(.*)_.*_counts.rds", "\\1", basename(counts_file))]

dat[, {
  sampleTable <- data.frame(unique(.SD[, .(cdition, 
                                           rep, 
                                           counts_file= basename(counts_file))]), 
                            row.names = "counts_file")
  sampleTable <- apply(sampleTable, 2, as.factor)
  # DESEq2 object
  if(!file.exists(dds_file)){
    DF <- .SD[, {
      .c <- data.table(readRDS(counts_file)$counts, keep.rownames = T)
      colnames(.c)[2] <- "counts"
      .c
    }, counts_file]
    DF <- data.frame(dcast(DF, rn~basename(counts_file), value.var = "counts"), row.names = "rn")
    DF <- DF[rowSums(DF)>10, match(rownames(sampleTable), colnames(DF))]
    dds <- DESeqDataSetFromMatrix(countData= DF, 
                                  colData= sampleTable, 
                                  design= ~rep+cdition)
    dds <- DESeq(dds)
    saveRDS(dds, dds_file)
  }else{
    dds <- readRDS(dds_file)
  }
  print(paste(dds_file, "DONE!"))
}, .(project, dds_file)]

#----------------------------------------------------------#
# Differential expression
#----------------------------------------------------------#
# Differential expression
diff <- unique(dat[, .(project, dds_file, cdition)])
diff <- diff[, CJ(nominator= cdition, denominator= cdition), .(dds_file, project)]
sel <- data.table(nominator= c("RNA_120hED",
                               "RNA_96hED",
                               "RNA_72hED",
                               "Ez18",
                               "Ez29",
                               "EzJ11",
                               "EzJ9",
                               "PH18",
                               "PH29",
                               "PHJ11",
                               "PHJ9",
                               "RNA_EZ7312A_ED",
                               "RNA_PCXT1092A_ED",
                               "RNA_PSCSUZ21B842D_ED",
                               "RNA_SUZ1212A_ED",
                               "RNA_PHRNAI_ED",
                               "RNA_4WED",
                               "RNA_8WED",
                               "RNA_14WED",
                               "RNA_ecdED"),
                  denominator= c("RNA_96hED",
                                 "RNA_72hED",
                                 "RNA_WTE1416",
                                 "W18",
                                 "W29",
                                 "WKD",
                                 "WKD",
                                 "W18",
                                 "W29",
                                 "WKD",
                                 "WKD",
                                 "RNA_2A_ED",
                                 "RNA_2A_ED",
                                 "RNA_42D_ED",
                                 "RNA_2A_ED",
                                 "RNA_WRNAI_ED",
                                 "RNA_WTED",
                                 "RNA_WTED",
                                 "RNA_WTED",
                                 "RNA_WTED"))
diff <- diff[sel, , on= c("nominator", "denominator")]
diff[, FC_file:= paste0("db/FC_tables/", project, "_", nominator, "_vs_", denominator, "_FC.txt"), .(nominator, denominator)]
diff[, {
  if(!file.exists(FC_file)){
    dds <- readRDS(dds_file)
    .c <- as.data.frame(lfcShrink(dds,
                                  type= "ashr",
                                  contrast= c("cdition", nominator, denominator)))
    # .c <- as.data.frame(results(dds,
    #                             contrast= c("cdition", V1, V2)))
    fwrite(.c, FC_file, col.names = T, row.names = T, sep ="\t", quote= F)
  }
  print(paste(FC_file, "DONE!"))
}, .(dds_file, FC_file, nominator, denominator)]

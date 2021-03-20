#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dds <- data.table(file= list.files("db/counts/", "_counts.rds$", recursive = T, full.names = T))
dds[, project:= tstrsplit(file, "/", keep= 4)]
dds[, dds_file:= paste0("db/dds/", project, "_dds.rds")]
dds[, rep:= as.numeric(sapply(gsub("_counts.rds", "", .SD$file), function(x) substr(x, nchar(x), nchar(x))))]
.expr <- "_1_counts.rds|_2_counts.rds|_3_counts.rds|_4_counts.rds|_rep1_counts.rds|_rep2_counts.rds|_rep3_counts.rds|_rep4_counts.rds"
dds[, cdition:= gsub(.expr, "", basename(file))]
dds[, cdition:= gsub("_$", "", cdition)]
dds[, {
  sampleTable <- data.frame(unique(.SD[, .(cdition, rep, counts_file= basename(file))]), row.names = "counts_file")
  sampleTable <- apply(sampleTable, 2, as.factor)
  if(!file.exists(dds_file)){
    # DESEq2
    DF <- .SD[, {
      .c <- data.table(readRDS(file)$counts, keep.rownames = T)
      colnames(.c)[2] <- "counts"
      .c
    }, file]
    DF <- data.frame(dcast(DF, rn~basename(file), value.var = "counts"), row.names = "rn")
    DF <- DF[rowSums(DF)>10, ]
    dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~rep+cdition)
    dds <- DESeq(dds)
    saveRDS(dds, dds_file)
  }else{
    dds <- readRDS(dds_file)
  }
  print(paste(dds_file, "DONE!"))
  # Differential expression
  diff <- CJ(sampleTable[, "cdition"], sampleTable[, "cdition"], unique = T)
  diff <- diff[V1!=V2]
  diff[, FC_file:= paste0("db/FC_tables/", project, "_", V1, "_vs_", V2, "_FC.txt"), .(V1, V2)]
  diff[, {
    if(!file.exists(FC_file)){
      .c <- as.data.frame(lfcShrink(dds, type= "ashr", contrast= c("cdition", V1, V2)))
      # .c <- as.data.frame(results(dds, contrast= c("cdition", V1, V2)))
      fwrite(.c, FC_file, col.names = T, row.names = T, sep ="\t", quote= F)
    }
    print(paste(FC_file, "DONE!"))
  }, .(FC_file, V1, V2)]
}, .(project, dds_file)]
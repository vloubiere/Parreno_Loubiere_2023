#----------------------------------------------------------#
# Compute counts
#----------------------------------------------------------#
counts <- data.table(file= list.files("db/bam/", ".bam$", recursive = T, full.names = T))
counts[, project:= tstrsplit(file, "/", keep= 4)]
counts[, {
  .c <- paste0("db/counts/", project)
  if(!dir.exists(.c)){
    dir.create(.c)
  }}, project]
counts[, counts_file:= paste0("db/counts/", project, "/", gsub(".bam", "_counts.rds", basename(file)))]
counts[, {
  if(!file.exists(counts_file)){
    if(.N==2){
      counts <- featureCounts(file, 
                              annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", 
                              isGTFAnnotationFile = T, 
                              isPairedEnd = T, 
                              nthreads = 8)
    }else if(.N==1){
      counts <- featureCounts(file, 
                              annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", 
                              isGTFAnnotationFile = T, 
                              isPairedEnd = F, 
                              nthreads = 8)
    }
    saveRDS(counts, counts_file)
  }
  print(paste(counts_file, "DONE!"))
}, .(file, counts_file)]


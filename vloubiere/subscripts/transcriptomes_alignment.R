#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
aln <- data.table(file= list.files("db/fastq", ".fq.gz", full.names = T, recursive = T))
aln[, project:= tstrsplit(file, "/", keep= 3)]
aln[, {
  .c <- paste0("db/bam/", project)
  if(!dir.exists(.c)){
    dir.create(.c)
  }}, project]
aln[, bam:= paste0("db/bam/", project, "/", .BY[1], ".bam"), gsub("_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
aln[, {
  if(!file.exists(bam)){
    print("START...")
    if(.N==2){
      stats <- capture.output(subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", 
                                      readfile1= file[1], 
                                      readfile2= file[2], 
                                      maxMismatches = 6, 
                                      nthreads = 10, 
                                      unique = T, 
                                      output_file= bam))
    }else if(.N==1){
      stats <- capture.output(subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", 
                                      readfile1= file[1], 
                                      nthreads = 10, 
                                      unique = T, 
                                      output_file= bam))
    }
    writeLines(stats, con = gsub(".bam$", "_stats.txt", bam))
  }
  print(paste(bam, "DONE!"))
}, bam]

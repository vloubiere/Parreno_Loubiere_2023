dat <- fread("external_data/filereport_read_run_PRJEB25187_tsv_tano.txt")
dat <- dat[grepl("^LGL_", sample_title)]

library(doParallel)
library(foreach)
Ncores <- detectCores()
cl <- makeCluster(Ncores)
registerDoParallel(cl)

registerDoParallel(cl) 
foreach(file= unlist(tstrsplit(dat$submitted_ftp, ";"))) %dopar% {
  download.file(paste0("ftp://", file), 
                destfile = paste0("db/fastq/DNA_tano/", basename(file)))
}

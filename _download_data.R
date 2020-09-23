setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions", ".R$", full.names = T), source)
require(data.table)
require(SRAdb)
require(GEOquery)
require(dplyr)

con <- dbConnect(SQLite(), "/groups/stark/vloubiere/projects/available_data/Rdata/SRAmetadb.sqlite")

#--------------------------------------------------------------#
# 1- Developmental RNA-Seq
#--------------------------------------------------------------#
# metadata <- rbindlist(lapply("GSE94915", getSRA, sra_con = con))
# metadata[, sample:= tstrsplit(experiment_title, " |;|_", keep= 3)]
# metadata[, rep:= tstrsplit(experiment_title, " |;|_", keep= 4)]
# metadata[rep=="repA", rep:= "rep1"]
# metadata[rep=="repB", rep:= "rep2"]
# metadata[rep=="repC", rep:= "rep3"]
# metadata <- metadata[, .(ftp= as.character(getFASTQinfo(in_acc = run, srcType = "ftp", sra_con = con)$ftp)), (metadata)]
# metadata[, filename:= paste0("db/fastq/dev_RNA/RNA_", sample, "_ED_", rep, ".fq.gz")]
# metadata[, bsub(paste("wget", ftp, "-O", filename)), (metadata)]

#--------------------------------------------------------------#
# 2- Ey PolII ChIP-Seq
#--------------------------------------------------------------#
# metadata <- rbindlist(lapply("GSE112868", getSRA, sra_con = con))
# metadata[, sample:= tstrsplit(experiment_title, " |;|_", keep= 2)]
# metadata <- metadata[sample %in% c("Ey", "No", "PolII")]
# metadata[, sample:= paste0(c("EY", "EY", "INPUT", "INPUT", "POLII"), "GSE112868")]
# metadata[, replicate:= tstrsplit(experiment_title, "Replicate |;", keep= 2)]
# metadata[5, replicate:= 1]
# metadata[, replicate:= paste0("rep", replicate)]
# metadata <- metadata[, .(ftp= as.character(getFASTQinfo(in_acc = run, srcType = "ftp", sra_con = con)$ftp)), (metadata)]
# metadata[, filename:= paste0("db/fastq/available_ChIP/ChIP_", sample, "_ED_", replicate, ".fq.gz")]
# metadata[, bsub(paste("wget", ftp, "-O", filename)), (metadata)]

#--------------------------------------------------------------#
# 3- ATAC-Seq
#--------------------------------------------------------------#
# metadata <- rbindlist(lapply("GSE59078", getSRA, sra_con = con))
# metadata <- metadata[grep("FRT82_ATAC-Seq", experiment_title)]
# metadata[, filename:= "ATAC_FRT82_ED_rep1.fq.gz"]
# metadata <- metadata[, .(ftp= as.character(getFASTQinfo(in_acc = run, srcType = "ftp", sra_con = con)$ftp)), (metadata)]
# metadata[, filename:= paste0("db/fastq/available_ATAC/", filename)]
# metadata[, bsub(paste("wget", ftp, "-O", filename)), (metadata)]






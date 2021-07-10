dat <- data.table(fq_file= c(list.files("F:/epiCancer/", recursive = T, pattern = ".fq.gz|.fastq.gz"),
                             list.files("db/fastq/", recursive = T, pattern = ".fq.gz|.fastq.gz")))
dat <- unique(dat[!grepl("merged", fq_file)])
dat[, project:= tstrsplit(fq_file, "/", keep= 1)]
dat[, type:= tstrsplit(project, "_", keep= 1)]
dat[, layout:= ifelse(.N==2, "PAIRED", "SINGLE"), gsub("_1.fq.gz|_2.fq.gz|_1.fastq.gz|_2.fastq.gz", "", fq_file)]

# DNA epicancer
dat[project=="DNA_epiCancer", cdition:= tstrsplit(basename(fq_file), "_H", keep = 1)]
dat[project=="DNA_epiCancer", rep:= as.character(1)]

# DNA Tano
meta_Tano <- fread("external_data/filereport_read_run_PRJEB25187_tsv_tano.txt")
cdition_Tano <- meta_Tano[dat[project=="DNA_tano", grep(patt, meta_Tano$submitted_ftp), .(patt= basename(fq_file))]$V1, sample_title]
dat[project=="DNA_tano", cdition:= cdition_Tano]
dat[project=="DNA_tano", rep:= 1]

# RNA epicancer
dat[project=="RNA_epiCancer", cdition:= gsub("_1.fq.gz|_2.fq.gz", "", basename(fq_file))]
dat[project=="RNA_epiCancer", rep:= gsub(".*_(.*)$", "\\1", cdition)]
dat[project=="RNA_epiCancer" & rep=="FACS", rep:= 1]
dat[project=="RNA_epiCancer", cdition:= gsub(paste0("_", rep, "$"), "", cdition), .(rep, cdition)]

# RNA development
dat[project=="RNA_development" & layout=="PAIRED", c("cdition", "rep"):= tstrsplit(gsub("_[1-9]{1}.fq.gz", "", basename(fq_file)), "_rep")]
dat[project=="RNA_development" & layout=="SINGLE", c("cdition", "rep"):= tstrsplit(gsub(".fq.gz", "", basename(fq_file)), "_rep")]

# RNA_mutants_SA2020
dat[project=="RNA_mutants_SA2020", c("cdition", "rep"):= tstrsplit(gsub(".fq.gz", "", basename(fq_file)), "_rep")]

# RNA_Paro_2018  
dat[project=="RNA_Paro_2018" & layout=="PAIRED", c("cdition", "rep"):= tstrsplit(gsub("_[1-9]{1}.fq.gz", "", basename(fq_file)), "_rep")]
dat[project=="RNA_Paro_2018" & layout=="SINGLE", c("cdition", "rep"):= tstrsplit(gsub(".fq.gz", "", basename(fq_file)), "_rep")]

# RNA_phRNAi_SA2020
dat[project=="RNA_phRNAi_SA2020", c("cdition", "rep"):= tstrsplit(gsub("_[1-9]{1}.fq.gz", "", basename(fq_file)), "_rep")]

dat[, prefix_output:= gsub("_1.fq.gz|_2.fq.gz|.fq.gz","",basename(fq_file))]
fwrite(dat, 
       "Rdata/raw_metadata.txt", 
       col.names = T, 
       row.names = F, 
       sep= "\t", 
       quote= F, 
       na= NA)

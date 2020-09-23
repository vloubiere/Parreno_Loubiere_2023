setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions", ".R$", full.names = T), source)
require(data.table)
require(SRAdb)
require(GEOquery)
require(dplyr)

#---------------------------------------------------#
# 1- Format metadata
#---------------------------------------------------#
# fq.gz
metadata <- data.table(file= list.files("/groups/stark/vloubiere/projects/epigenetic_cancer/db/fastq", ".fq.gz", full.names = T, recursive = T))
# method and projects
metadata[, method:= tstrsplit(file, "/|_", keep= 11)]
metadata[, project:= tstrsplit(file, "/|_", keep= 10)]
# Sample, replicates and IDs
metadata[project=="SA2020", c("sample", "rep"):= tstrsplit(basename(file), "_|.fq", keep= c(2, 4))]
metadata[project=="available", c("sample", "rep"):= tstrsplit(basename(file), "_|.fq", keep= c(2, 4))]
metadata[project=="dev", c("sample", "rep"):= tstrsplit(basename(file), "_|.fq", keep= c(2, 3))]
metadata[project=="epiCancer", c("sample", "rep"):= tstrsplit(basename(file), "_|.fq", keep= c(1,2))]
metadata[, id:= paste0(sample, "_", rep), .(sample, rep)]
# Genomes and prefixes
metadata[method=="RNA" & (project=="epiCancer" | grepl("RNAI$", sample)), genome:= "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_150bp/"]
metadata[method=="RNA" & is.na(genome) & project=="dev" & grepl("hED$", sample), genome:= "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_50bp/"]
metadata[method=="RNA" & is.na(genome) & project=="dev" & grepl("E1416$", sample), genome:= "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_75bp/"]
metadata[method=="RNA" & is.na(genome) & project=="SA2020", genome:= "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_50bp/"]
metadata[method=="RNA", bam_prefix:= paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bam/", project, "_RNA/", id), project]
metadata[method=="RNA", counts_prefix:= paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/counts/", project, "/", id), project]
metadata[method=="RNA", bed_prefix:= paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bed/RNA/", id), project]
metadata[method=="RNA", bw_prefix:= paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/RNA/", id), project]

fwrite(metadata, "db/metadata/metadata_all.txt", col.names = T, row.names = F, sep= "\t", quote= F, na= NA)

#---------------------------------------------------#
# 2- commands
#---------------------------------------------------#
metadata <- metadata[, .(cmd= 
  if(method=="ChIP")
  {
  
  c("module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript /groups/stark/vloubiere/pipelines/ChIPseq_pipeline.R", 
    file[1], 
    ifelse(.N==2, file[2], NA), 
    "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome", 
    "dm6", 
    paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bed/ChIP/", id), 
    paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ChIP/", id))
  
  }else if(method=="ATAC"){
    
    c("module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript /groups/stark/vloubiere/pipelines/ChIPseq_pipeline.R", 
      file[1], 
      ifelse(.N==2, file[2], NA), 
      "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome", 
      "dm6", 
      paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bed/ATAC/", id), 
      paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/db/bw/ATAC/", id))
    
  }else if(method=="RNA"){
    
    c("module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript /groups/stark/vloubiere/pipelines/RNAseq_pipeline.R", 
      file[1], 
      ifelse(.N==2, file[2], NA), 
      genome,
      bam_prefix, 
      bed_prefix, 
      bw_prefix, 
      counts_prefix, 
      "/groups/stark/vloubiere/genomes/flybase/dmel-all-r6.35_simplified.gtf")
    
  }), .(id, method, project, genome, bam_prefix, bed_prefix, bw_prefix, counts_prefix)]
metadata <- metadata[, .(cmd= paste0(cmd, collapse= " ")), .(id, method, project)]

#---------------------------------------------------#
# 3- RUN
#---------------------------------------------------#
logs <- "/groups/stark/vloubiere/projects/epigenetic_cancer/db/logs/"
metadata[grepl("E1416", id), bsub(cmd, o= paste0(logs, id), e= paste0(logs, id), cores = 12, m = 40, name = paste0("vl_", id)), .(cmd, id)]
# metadata[project %in% c("available", "dev"), bsub(cmd, o= paste0(logs, id), e= paste0(logs, id), cores = 12, m = 40, name = paste0("vl_", id)), .(cmd, id)]
# metadata[, bsub(cmd, o= paste0(logs, id), e= paste0(logs, id), cores = 12, m = 40, name = paste0("vl_", id)), .(cmd, id)]

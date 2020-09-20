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
metadata <- data.table(file= list.files("/groups/stark/vloubiere/data/epigenetic_cancer", ".fq.gz", full.names = T, recursive = T))
# method and projects
metadata[, method:= if(grepl("_RNA$", dirname(file))){"RNA"}else if(grepl("_ChIP$", dirname(file))){"ChIP"}, file]
metadata[, project:= if(grepl("/SA2020_", dirname(file))){"SA2020"}else if(grepl("/epiCancer_", dirname(file))){"epiCancer"}, file]
# Sample, replicates and IDs
metadata[project=="SA2020", c("sample", "rep"):= tstrsplit(basename(file), "_|.fq", keep= c(2, 4))]
metadata[project=="epiCancer", c("sample", "rep"):= tstrsplit(basename(file), "_|.fq", keep= c(1,2))]
metadata[, id:= paste0(sample, "_", rep), .(sample, rep)]

#---------------------------------------------------#
# 2- commands
#---------------------------------------------------#
# ChIP
ChIP_pipeline <- "module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript /groups/stark/vloubiere/pipelines/ChIPseq_pipeline.R"
bwt2_idx <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/Bowtie2Index/genome"
bed_prefix <- "/groups/stark/vloubiere/data/epigenetic_cancer/bed/"
bw_prefix <- "/groups/stark/vloubiere/data/epigenetic_cancer/bw/"

RNA_pipeline <- "module load r/3.6.2-foss-2018b; /software/2020/software/r/3.6.2-foss-2018b/bin/Rscript /groups/stark/vloubiere/pipelines/RNAseq_pipeline.R"
genome150 <- "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_150bp/"
genome50 <- "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_50bp/"
bam_prefix_SA <- "/groups/stark/vloubiere/data/epigenetic_cancer/bam/SA2020_RNA/"
bam_prefix_Epi <- "/groups/stark/vloubiere/data/epigenetic_cancer/bam/epiCancer_RNA/"
counts_prefix_SA <- "/groups/stark/vloubiere/data/epigenetic_cancer/counts/SA2020/"
counts_prefix_Epi <- "/groups/stark/vloubiere/data/epigenetic_cancer/counts/epiCancer/"
gtf <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf"

metadata <- metadata[, .(cmd=
                           {
                             if(method=="ChIP")
                             {
                               c(ChIP_pipeline, file[1], ifelse(.N==2, file[2], NA), bwt2_idx, "dm6", paste0(bed_prefix, id), paste0(bw_prefix, id))
                             }
                             if(method=="RNA")
                             {
                               bam_prefix <- paste0(ifelse(project=="SA2020", bam_prefix_SA, bam_prefix_Epi), id)
                               counts_prefix <- paste0(ifelse(project=="SA2020", counts_prefix_SA, counts_prefix_Epi), id)
                               genome <- ifelse(.N==2, genome150, genome50)
                               c(RNA_pipeline, file[1], ifelse(.N==2, file[2], NA), genome, bam_prefix, counts_prefix, gtf)
                             }
                           }), .(id, method, project)]
metadata <- metadata[, .(cmd= paste0(cmd, collapse= " ")), .(id, method, project)]

#---------------------------------------------------#
# 3- RUN
#---------------------------------------------------#
logs <- "/groups/stark/vloubiere/data/epigenetic_cancer/logs/"
metadata[, bsub(cmd, o= logs, e= logs, cores = 8, m = 20), cmd]






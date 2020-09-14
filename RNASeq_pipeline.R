setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
source("/groups/stark/vloubiere/scripts/R_functions/bsub_wrapper_1.0.R")
require(data.table)

metadata <- fread("/groups/stark/vloubiere/data/epigenetic_cancer/metadata/metadata.txt")
metadata[, fq1:= list.files("/groups/stark/vloubiere/data/epigenetic_cancer/fastq", paste0("^", ID, "*_1.fq.gz"), full.names = T), ID]
metadata[, fq2:= list.files("/groups/stark/vloubiere/data/epigenetic_cancer/fastq", paste0("^", ID, "*_2.fq.gz"), full.names = T), ID]

# Generate sam file
metadata[, bam_file:= paste0("/groups/stark/vloubiere/data/epigenetic_cancer/bam/", ID, "Aligned.out.bam"), ID]
metadata[!file.exists(bam_file), star_cmd:= paste("module load star/2.7.1a-foss-2018b;
                                                  /software/2020/software/star/2.7.1a-foss-2018b/bin/STAR
                                                  --runThreadN 8
                                                  --genomeDir /groups/stark/vloubiere/genomes/STAR_genome_150bp/
                                                  --outSAMtype BAM Unsorted
                                                  --outSJfilterReads Unique
                                                  --outFilterMultimapNmax 1
                                                  --readFilesCommand zcat
                                                  --outFileNamePrefix", paste0("/groups/stark/vloubiere/data/epigenetic_cancer/bam/", ID),
                                                 "--readFilesIn", fq1, fq2), c(colnames(metadata))]
 
metadata[, counts_file:= paste0("/groups/stark/vloubiere/projects/epigenetic_cancer/counts/", ID, "_counts.txt"), bam_file]
metadata[!file.exists(counts_file),
         {
           htseq_cmd:= paste("module load  htseq/0.11.2-foss-2018b-python-3.6.6;
                             /software/2020/software/htseq/0.11.2-foss-2018b-python-3.6.6/bin/htseq-count -s no -f bam", bam_file, 
                            "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf > ", counts_file)
         }, .(bam_file, counts_file)]

RUN <- melt(metadata, id.vars = "ID", measure.vars = patterns("cmd$"))
RUN <- RUN[!is.na(value), .(cmd= paste(value, collapse = ";")), ID]
# RUN[, bsub(cmd, name = ID, cores = 8, m = 16), ID]
RUN[, bsub(cmd, name = ID), ID]
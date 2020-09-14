source("/groups/stark/vloubiere/scripts/R_functions/bsub_wrapper_1.0.R")

cores <- 6
genomeDir <- "/groups/stark/vloubiere/genomes/STAR_genome_150bp/"
dm3.fa <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa"
dm3_genes.gtf <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf"
cmd <- paste0("module load star/2.7.1a-foss-2018b; /software/2020/software/star/2.7.1a-foss-2018b/bin/STAR 
                --runMode genomeGenerate 
                --runThreadN ", cores,
              " --genomeDir ", genomeDir, 
              " --genomeFastaFiles ", dm3.fa,
              " --sjdbGTFfile ", dm3_genes.gtf , " --sjdbOverhang 149")
bsub(cmd, cores= cores)

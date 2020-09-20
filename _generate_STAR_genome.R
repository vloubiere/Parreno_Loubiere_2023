source("/groups/stark/vloubiere/functions/bsub_wrapper.R")

#-------------------------------------------------------------#
# dm3
#-------------------------------------------------------------#
# cores <- 6
# genomeDir <- "/groups/stark/vloubiere/genomes/STAR_genome/dm3/STAR_genome_150bp/"
# dm3.fa <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm3/Sequence/WholeGenomeFasta/genome.fa"
# dm3_genes.gtf <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm3/Annotation/Genes/genes.gtf"
# cmd <- paste0("module load star/2.7.1a-foss-2018b; /software/2020/software/star/2.7.1a-foss-2018b/bin/STAR 
#                 --runMode genomeGenerate 
#                 --runThreadN ", cores,
#               " --genomeDir ", genomeDir, 
#               " --genomeFastaFiles ", dm3.fa,
#               " --sjdbGTFfile ", dm3_genes.gtf , " --sjdbOverhang 149")
# bsub(cmd, cores= cores)

#-------------------------------------------------------------#
# dm6
#-------------------------------------------------------------#
# 150 nt
cores <- 6
genomeDir <- "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_150bp/"
dm6.fa <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"
dm6_genes.gtf <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf"
cmd <- paste0("module load star/2.7.1a-foss-2018b; /software/2020/software/star/2.7.1a-foss-2018b/bin/STAR 
                --runMode genomeGenerate 
                --runThreadN ", cores,
              " --genomeDir ", genomeDir, 
              " --genomeFastaFiles ", dm6.fa,
              " --sjdbGTFfile ", dm6_genes.gtf , " --sjdbOverhang 149")
bsub(cmd, o= genomeDir, e= genomeDir, cores= cores)
# 50nt
cores <- 6
genomeDir <- "/groups/stark/vloubiere/genomes/STAR_genome/dm6/STAR_genome_50bp/"
dm6.fa <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa"
dm6_genes.gtf <- "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Annotation/Genes/genes.gtf"
cmd <- paste0("module load star/2.7.1a-foss-2018b; /software/2020/software/star/2.7.1a-foss-2018b/bin/STAR 
                --runMode genomeGenerate 
                --runThreadN ", cores,
              " --genomeDir ", genomeDir, 
              " --genomeFastaFiles ", dm6.fa,
              " --sjdbGTFfile ", dm6_genes.gtf , " --sjdbOverhang 49")
bsub(cmd, o= genomeDir, e= genomeDir, cores= cores)

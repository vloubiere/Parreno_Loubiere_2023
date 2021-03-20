#----------------------------------------------------------#
# Build dm6 index
#----------------------------------------------------------#
setwd("/_R_data/genomes/dm6/subreadr_index/")
ref <- "/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
buildindex(basename= "subreadr_dm6_index", reference= ref)
setwd("/_R_data/projects/epigenetic_cancer/")
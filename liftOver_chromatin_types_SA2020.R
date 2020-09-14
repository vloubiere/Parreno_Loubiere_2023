setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
source("/groups/stark/vloubiere/scripts/R_functions/bedtools_wrappers_1.0.R")
require(circlize)
require(pheatmap)
require(data.table)
require(gridExtra)
require(clusterProfiler)
require(org.Dm.eg.db)
require(rtracklayer)
require(colorRamps)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)

# Chromatin type
# # Lift over dm6 ->3
# chrom <- import("/groups/stark/vloubiere/projects/SA_2019/data_SA_2019/chrom_types_ED_SA_2020_dm6.txt", format = "BED")
# chain <- import("/groups/stark/vloubiere/genomes/dm6ToDm3.over.chain")
# frag <- sample(20, length(chrom), replace = T)
# res <- list()
# for(i in seq(20))
# {
#   res[[i]] <- Reduce(c, liftOver(chrom[frag==i], chain))
# }
# res <- sortSeqlevels(Reduce(c, res))
# res <- sort(res)
# export(res, "chrom_types/chromatin_types_SA_2020_dm3.bed")

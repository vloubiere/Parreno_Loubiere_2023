setwd("/groups/stark/vloubiere/genomes/flybase/")
sapply(list.files("/groups/stark/vloubiere/functions", ".R$", full.names = T), source)

require(rtracklayer)

gtf <- import("/groups/stark/vloubiere/genomes/flybase/dmel-all-r6.35.gtf")
gtf <- gtf[as.character(seqnames(gtf)) %in% c("2L", "2R", "3L", "3R", "4", "X", "Y")]
seqlevelsStyle(gtf) <- "UCSC"
export(gtf, "/groups/stark/vloubiere/genomes/flybase/dmel-all-r6.35_simplified.gtf")

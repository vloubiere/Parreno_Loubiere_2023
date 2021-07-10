require(data.table)
require(gridExtra)

stats <- data.table(file= list.files("db/bam/", ".bam.summary", full.names = T, recursive = T))
stats[, cdition:= gsub(".bam.summary", "" , basename(file))]
stats <- stats[, fread(file), (stats)]
stats <- dcast(stats, cdition~V1, value.var = "V2")
stats[is.na(Total_reads), Total_reads:= Total_fragments]
stats[is.na(Uniquely_mapped_reads), Uniquely_mapped_reads:= Uniquely_mapped_fragments]
stats[, Percentage:= Uniquely_mapped_reads/Total_reads*100, .(Total_reads, Uniquely_mapped_reads)]
stats[, Uniquely_mapped_reads:= formatC(Uniquely_mapped_reads, big.mark = ",", format = "d")]

pdf("pdf/Alignment_statistics.pdf", height = 25)
grid.table(stats[, .(cdition, Uniquely_mapped_reads, Percentage)])
dev.off()

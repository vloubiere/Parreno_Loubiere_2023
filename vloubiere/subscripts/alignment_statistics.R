stats <- data.table(file= list.files("db/bam/", ".bam.summary", full.names = T, recursive = T))
stats[, cdition:= gsub(".bam.summary", "" , basename(file))]
stats <- stats[, fread(file), (stats)]
stats <- dcast(stats, cdition~V1, value.var = "V2")
stats[is.na(Total_reads), Total_reads:= Total_fragments]
stats[is.na(Uniquely_mapped_reads), Uniquely_mapped_reads:= Uniquely_mapped_fragments]
stats[, Percentage:= Uniquely_mapped_reads/Total_reads*100, .(Total_reads, Uniquely_mapped_reads)]
stats[, Uniquely_mapped_reads:= formatC(Uniquely_mapped_reads, big.mark = ",", format = "d")]

# stats <- data.table(file= list.files("db/bam/", "stats.txt", full.names = T, recursive = T))
# stats <- stats[, .(Mapped_reads= readLines(file)), file]
# stats <- stats[grepl("Mapped :", Mapped_reads)]
# stats[, Mapped_reads:= gsub("\\||Mapped| |:|fragments|, wherein", "", Mapped_reads)]
# stats[, Mapped_reads:= gsub("reads", "", Mapped_reads)]
# stats[, c("Mapped_reads", "Percentage"):= tstrsplit(Mapped_reads, "\\(|\\)", keep= c(1,2))]
# stats[, Mapped_reads:= formatC(as.numeric(Mapped_reads), big.mark = ",", format = "d")]
# stats[, cdition:= gsub("_stats.txt", "" , basename(file))]

pdf("pdf/Alignment_statistics.pdf", height = 25)
grid.table(stats[, .(cdition, Uniquely_mapped_reads, Percentage)])
dev.off()

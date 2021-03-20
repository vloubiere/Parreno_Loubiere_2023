stats <- data.table(file= list.files("db/bam/", "stats.txt", full.names = T))
stats <- stats[grepl("^Ez|^PH|^W", basename(file))]
stats <- stats[, .(Mapped_reads= readLines(file)), file]
stats <- stats[grepl("Mapped :", Mapped_reads)]
stats[, Mapped_reads:= gsub("\\||Mapped| |:|fragments|, wherein", "", Mapped_reads)]
stats[, Mapped_reads:= gsub("\\(", " (", Mapped_reads)]
stats[, file:= gsub("_stats.txt", "" , basename(file))]

pdf("pdf/Alignment_statistics.pdf", height = 10)
grid.table(stats)
dev.off()

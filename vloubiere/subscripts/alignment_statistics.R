require(data.table)
require(gridExtra)

stats <- data.table(file= list.files(c("db/bam/", 
                                       "/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/"), 
                                     ".bam.summary", 
                                     full.names = T, 
                                     recursive = T))
stats[, cdition:= tstrsplit(file, "/bam//|.bam.summary", keep = 2)]
stats <- stats[, fread(file), (stats)]
stats <- dcast(stats, cdition~V1, value.var = "V2")
stats[is.na(Total_reads), Total_reads:= Total_fragments]
stats[is.na(Uniquely_mapped_reads), Uniquely_mapped_reads:= Uniquely_mapped_fragments]
stats[, Percentage:= Uniquely_mapped_reads/Total_reads*100, .(Total_reads, Uniquely_mapped_reads)]
stats[, Uniquely_mapped_reads:= formatC(Uniquely_mapped_reads, big.mark = ",", format = "d")]

pdf("pdf/alignment/Alignment_statistics.pdf", height = 35, width = 12)
grid.table(stats[, .(cdition, Uniquely_mapped_reads, Percentage)])
dev.off()

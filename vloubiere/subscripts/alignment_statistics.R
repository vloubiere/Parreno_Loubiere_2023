require(data.table)
require(gridExtra)

stats <- fread("Rdata/processed_metadata_RNA.txt")
stats <- stats[grepl("^epiCancer", DESeq2_object), fread(gsub(".bam$", ".bam.summary", bam)), .(bam, cdition, rep, DESeq2_object)]
stats <- dcast(stats[V1 %in% c("Total_fragments", "Uniquely_mapped_fragments")], 
               cdition+rep+DESeq2_object~V1, 
               value.var = "V2")
stats[, Percentage:= round(Uniquely_mapped_fragments/Total_fragments*100, 1)]
cols <- c("Total_fragments", "Uniquely_mapped_fragments")
stats[, (cols):= lapply(.SD, formatC, big.mark = ",", format = "d"), .SDcols= cols]
setkeyv(stats, c("DESeq2_object", "cdition", "rep"))
setcolorder(stats, "DESeq2_object")

pdf("pdf/Figures/Alignment_statistics_RNA.pdf", 
    height = 35, 
    width = 12)
grid.table(stats)
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(gridExtra)

stats <- fread("Rdata/processed_metadata_CUTNRUN.txt")
stats <- stats[, fread(gsub(".bam$", ".bam.summary", bam)), .(bam, cdition, rep, system)]
stats <- dcast(stats[V1 %in% c("Total_fragments", "Uniquely_mapped_fragments")], 
               cdition+rep+system~V1, 
               value.var = "V2")
stats[, Percentage:= round(Uniquely_mapped_fragments/Total_fragments*100, 1)]
cols <- c("Total_fragments", "Uniquely_mapped_fragments")
stats[, (cols):= lapply(.SD, formatC, big.mark = ",", format = "d"), .SDcols= cols]
setkeyv(stats, c("system", "cdition", "rep"))
setcolorder(stats, "system")

pdf("pdf/Figures/Alignment_statistics_RNA.pdf", 
    height = 12, 
    width = 10)
grid.table(stats)
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(gridExtra)

stats <- fread("Rdata/processed_metadata_CUTNRUN.txt")
stats[, summary:= paste0(bam, ".summary.txt")]
# Compute stats
mclapply(stats[!file.exists(summary), bam], function(x) 
  system(paste0("samtools stats ", x, " | grep ^SN > ", x, ".summary.txt")), 
  mc.preschedule = T)
stats[, Total_fragments:= fread(cmd= paste0("zcat ", fq1, " | wc -l "))[, V1]/4, fq1]
stats[, Uniquely_mapped_fragments:= fread(paste0(bam, ".summary.txt"))[7, V3]/2, bam]
stats[, Percentage:= round(Uniquely_mapped_fragments/Total_fragments*100, 1)]
cols <- c("Total_fragments", "Uniquely_mapped_fragments")
stats[, (cols):= lapply(.SD, formatC, big.mark = ",", format = "d"), .SDcols= cols]

pdf("pdf/cutnrun_alignment_statistics_cutnrun.pdf", 
    height = 20, 
    width = 10)
grid.table(stats[, .(ChIP, cdition, rep, Total_fragments, Uniquely_mapped_fragments, Percentage)])
dev.off()
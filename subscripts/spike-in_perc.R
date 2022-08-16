setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)

# SPIKE IN fraction
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
meta[, ChIP_reads:= fread(cmd= paste0("wc -l ", ChIP_bed))$V1, ChIP_bed]
meta[, spikein_reads:= fread(cmd= paste0("wc -l ", spikein_bed))$V1, spikein_bed]
pl <- unique(meta[, .(name= paste0(ChIP, "_", cdition, "_", rep), ChIP_reads, spikein_reads)])
pl <- t(as.matrix(pl, 1))

pdf("pdf/cutnrun/spikein_percentage.pdf", width = 10)
layout(matrix(c(1,3,2), nrow= 3), 
       heights = c(1.7,0.2,2))
Cc <- c("grey20", "grey80")
par(mar= c(2,5,2,2))
barplot(pl, 
        beside = T,
        las= 2, 
        col= Cc,
        ylab= "N reads\n", 
        xaxt= "n")
legend("topright",
       fill= Cc,
       legend= c("ChIP", "Spike-in"),
       bty= "n",
       xpd= T)
par(mar= c(12,5,2,2))
bar <- barplot(pl["spikein_reads",]/pl["ChIP_reads",]*100, 
               ylim= c(0, 2.5),
               beside = T,
               las= 2,
               ylab= "% spike-in reads")
text(bar, 
     grconvertY(1, "nfc", "user"),
     paste0(round(pl["spikein_reads",]/pl["ChIP_reads",]*100, 2), "%"),
     pos= 1,
     xpd= T)
dev.off()
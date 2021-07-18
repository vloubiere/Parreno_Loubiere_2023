setwd("D:/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(dplyr)
require(vlfunctions)

#--------------------------------------------#
# ATAC-Seq peak calling
#--------------------------------------------#
ATAC <- list.files("D:/_R_data/projects/public_data/dm6/bam/",
                   paste0(dat[grepl("ATAC", type), Run], collapse= "|"),
                   full.names = T)
ATAC <- ATAC[grepl(".bed$", ATAC)]
ATAC_peaks_ED <- vl_peakCalling(ChIP_bed = ATAC, 
                                BSgenome = BSgenome.Dmelanogaster.UCSC.dm6, 
                                min_N_replicates = 8)

dir.create("db/peaks", 
           showWarnings = F)
fwrite(ATAC_peaks_ED, 
       "db/peaks/ATAC_seq_peaks_ED.txt", 
       col.names = T, 
       row.names = F, 
       sep= "\t", 
       quote=F)

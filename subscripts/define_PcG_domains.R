setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# List bw files of interest ----
dat <- data.table(bw= c("db/bw/cutnrun/H3K27me3_PH18_merge.bw",
                        "db/bw/cutnrun/H3K27me3_PH29_merge.bw",
                        "db/bw/cutnrun/H3K27me3_PHD11_merge.bw"),
                  domains= c("db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak",
                             "db/peaks/cutnrun/H3K27me3_PH29_confident_peaks.broadPeak",
                             "db/peaks/cutnrun/H3K27me3_PHD11_confident_peaks.broadPeak"),
                  cdition= c("PH18", "PH29", "PHD11"))
dat[, merged_domains:= paste0("db/bed/merged_K27_domains/", cdition, ".bed")]

# Bin the genome to find low complexity regions with no signal
bin <- vl_binBSgenome(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6,
                      bins.width = 500,
                      steps.width = 250,
                      restrict.seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"))
.q <- vl_bw_coverage(bin, bw = "db/bw/cutnrun/H3K27me3_PH18_merge.bw")
gaps <- vl_collapseBed(bin[.q==0])

# Compute merged K27me3 domains ----
dat[, {
  # Import K27 broadpeaks and merge with gaps
  cmb <- rbindlist(list(gap= gaps,
                        K27= vl_importBed(domains)),
                   fill = T,
                   idcol = "class")
  cmb$idx <- vl_collapseBed(cmb,
                            ming.gap = 2.5e3,
                            return.idx.only = T)
  cmb <- cmb[idx %in% cmb[class=="K27", idx]]
  merge <- vl_collapseBed(cmb,
                          ming.gap = 2.5e3)
  # Save 
  fwrite(merge,
         merged_domains,
         col.names = F,
         row.names = F,
         sep= "\t",
         quote= F)
  print(paste(.GRP, .NGRP, sep= " / "))
}, .(cdition, domains, merged_domains)]

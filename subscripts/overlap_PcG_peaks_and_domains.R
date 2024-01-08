setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import PH18 domains and peaks
ph.18 <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")[signalValue>3]
k27.18 <- vl_importBed("db/bed/merged_K27_domains/PH18.bed")
k118.18 <- vl_importBed("db/peaks/cutnrun/H2AK118Ub_PH18_confident_peaks.broadPeak")
k27.18 <- vl_intersectBed(k27.18, ph.18)
k118.18 <- vl_intersectBed(k118.18, ph.18)
ph.18 <- vl_intersectBed(ph.18, rbind(k27.18[, .(seqnames, start, end)], k118.18[, .(seqnames, start, end)]))

# Intersect per condition
k27.29 <- vl_intersectBed(k27.18,
                          "db/bed/merged_K27_domains/PH29.bed",
                          min.overlap = k27.18[, end-start+1]/2)
k27.d11 <- vl_intersectBed(k27.18,
                           "db/bed/merged_K27_domains/PHD11.bed",
                           min.overlap = k27.18[, end-start+1]/2)
ph.29 <- vl_intersectBed(ph.18,
                         "db/peaks/cutnrun/PH_PH29_confident_peaks.narrowPeak")
ph.d11 <- vl_intersectBed(ph.18,
                          "db/peaks/cutnrun/PH_PHD11_confident_peaks.narrowPeak")
k118.29 <- vl_intersectBed(k118.18,
                           "db/peaks/cutnrun/H2AK118Ub_PH29_confident_peaks.broadPeak",
                           min.overlap = k118.18[, end-start+1]/2)
k118.d11 <- vl_intersectBed(k118.18,
                            "db/peaks/cutnrun/H2AK118Ub_PHD11_confident_peaks.broadPeak",
                            min.overlap = k118.18[, end-start+1]/2)

# Plots
pdf("pdf/review_overlap_PcG_features.pdf", 7, 2.75)
vl_par(mfrow= c(1,3),
       mai= c(1, 1.5, .7, .2))
vl_upset_plot(list("No ph-KD"= ph.18$name,
                   "Constant ph-KD"= ph.29$name,
                   "Transient ph-KD"= ph.d11$name),
              cex.grid = 1.2,
              grid.hex = 1.2)
title(main= "PH")
vl_upset_plot(list("No ph-KD"= k27.18$name,
                   "Constant ph-KD"= k27.29$name,
                   "Transient ph-KD"= k27.d11$name),
              cex.grid = 1.2,
              grid.hex = 1.2)
title(main= "H3K27me3")
vl_upset_plot(list("No ph-KD"= k118.18$name,
                   "Constant ph-KD"= k118.29$name,
                   "Transient ph-KD"= k118.d11$name),
              cex.grid = 1.2,
              grid.hex = 1.2,
              show.empty = T)
title(main= "H2AK118Ub")
dev.off()
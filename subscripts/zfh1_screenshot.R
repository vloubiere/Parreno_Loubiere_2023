require(vlfunctions)

sel <- as.data.table(GRanges("chr3R:30745267-30805712"))

cutoffPH <- 10
tmp1 <- tempfile(fileext = ".narrowPeak")
fwrite(vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")[signalValue>cutoffPH],
       tmp1, 
       col.names = F, 
       row.names = F,
       sep= "\t")


tracks <- c("db/bw/cutnrun/PH_PH18_merge.bw", 
            tmp1,
            "db/bw/cutnrun/H2AK118Ub_PH18_merge.bw",
            "db/peaks/cutnrun/H2AK118Ub_PH18_confident_peaks.broadPeak",
            "db/bw/cutnrun/H3K27me3_PH18_merge.bw",
            "db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak",
            "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PH18_confident_peaks.narrowPeak")

pdf("pdf/Extended_data_5_screenshot_zfh1.pdf",
    width = 4, 
    height = 2.5)
par(mar= c(6.1, 6.1, 1.1, 2.1))
vl_screenshot(bed = sel,
              tracks= tracks, 
              names = c("PH", "", "H2AK118Ub", "", "H3K27me3", "", "H3K27Ac", ""), 
              max= c(100, 10, 20, 10),
              col= rep(c("grey0", "grey20", "grey40", "grey60"), each= 2), 
              widths = c(150L, 20L),
              genome = "dm6")
dev.off()

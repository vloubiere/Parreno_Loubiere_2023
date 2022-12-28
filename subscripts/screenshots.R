require(vlfunctions)

sel <- c("chrX:18230628-18329654",
         "chr3R:30735267-30797712",
         "chr2R:11450419-11590759",
         "chr2L:7275066-7335276")
sel <- as.data.table(GRanges(sel))

cutoffPH <- 10
tmp1 <- tempfile(fileext = ".narrowPeak")
fwrite(vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")[signalValue>cutoffPH],
       tmp1, 
       col.names = F, 
       row.names = F,
       sep= "\t")
tmp2 <- tempfile(fileext = ".narrowPeak")
fwrite(vl_importBed("db/peaks/cutnrun/PH_PHD11_confident_peaks.narrowPeak")[signalValue>cutoffPH],
       tmp2, 
       col.names = F, 
       row.names = F,
       sep= "\t")


tracks <- c("db/bw/cutnrun/PH_PH18_merge.bw", 
            tmp1,
            "db/bw/cutnrun/PH_PH29_merge.bw",
            "db/peaks/cutnrun/PH_PH29_confident_peaks.narrowPeak",
            "db/bw/cutnrun/PH_PHD11_merge.bw",
            tmp2,
            "db/bw/cutnrun/H3K27me3_PH18_merge.bw",
            "db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak",
            "db/bw/cutnrun/H3K27me3_PH29_merge.bw",
            "db/peaks/cutnrun/H3K27me3_PH29_confident_peaks.broadPeak",
            "db/bw/cutnrun/H3K27me3_PHD11_merge.bw",
            "db/peaks/cutnrun/H3K27me3_PHD11_confident_peaks.broadPeak",
            "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PH18_confident_peaks.narrowPeak",
            "db/bw/cutnrun/H3K27Ac_PH29_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PH29_confident_peaks.narrowPeak",
            "db/bw/cutnrun/H3K27Ac_PHD11_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PHD11_confident_peaks.narrowPeak")
names <- c("No ph-KD", "", "Constant ph-KD", "", "Transient ph-KD", "")

pdf("pdf/Figure_3_screenshot_reversible_irreversible_genes.pdf",
    width = 10, 
    height = 4.5)
par(mar= c(6.1, 6.1, 1.1, 2.1))
vl_screenshot(bed = sel,
              tracks= tracks, 
              names = names, 
              max= c(100, 100, 350, 22, 22, 10, 10, 10, 10),
              col= rep(c("grey0", "grey35", "grey70"), each= 6), 
              widths = c(150L, 20L),
              genome = "dm6")
dev.off()

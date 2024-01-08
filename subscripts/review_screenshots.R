setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)

sel <- data.table(seqnames= c("chrX",
                              "chr3R",
                              "chr2R",
                              "chr2L"),
                  start= c(18230628,
                           30735267,
                           11450419,
                           7275066),
                  end= c(18329654,
                         30797712,
                         11590759,
                         7335276))
sel[, c("start", "end"):= .(start-((end-start)/5), end+((end-start)/5))]

# Select strong PH peaks ----
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

# bw tracks ---
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
            "db/bw/cutnrun/H2AK118Ub_PH18_merge.bw",
            "db/peaks/cutnrun/H2AK118Ub_PH18_confident_peaks.broadPeak",
            "db/bw/cutnrun/H2AK118Ub_PH29_merge.bw",
            "db/peaks/cutnrun/H2AK118Ub_PH29_confident_peaks.broadPeak",
            "db/bw/cutnrun/H2AK118Ub_PHD11_merge.bw",
            "db/peaks/cutnrun/H2AK118Ub_PHD11_confident_peaks.broadPeak",
            "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PH18_confident_peaks.narrowPeak",
            "db/bw/cutnrun/H3K27Ac_PH29_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PH29_confident_peaks.narrowPeak",
            "db/bw/cutnrun/H3K27Ac_PHD11_merge.bw",
            "db/peaks/cutnrun/H3K27Ac_PHD11_confident_peaks.narrowPeak")
names <- c("No ph-KD", "", "Constant ph-KD", "", "Transient ph-KD", "")

# Plot ----
pdf("pdf/review_screenshot_reversible_irreversible_genes.pdf",
    width = 10, 
    height = 4.5)
par(mar= c(6.1, 6.1, 1.1, 2.1),
    mgp= c(-1, -1, 0),
    cex= 1,
    las= 1,
    cex.axis= 6/12,
    cex.lab= 9/12)
vl_screenshot(bed = sel,
              tracks= tracks, 
              names = names, 
              max= c(100, 100, 350, 18, 24, 15, 15, 15, 15, 10, 10, 10),
              col= rep(c("grey0", "grey40", "grey0", "grey40"), each= 6), 
              widths = c(150L, 20L),
              genome = "dm6")
dev.off()
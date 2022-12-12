setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# Import genes
genes <- fread("Rdata/final_gene_features_table.txt")
proms <- vl_resizeBed(genes, "start", 0, 0)

# PH peaks
PH <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
setorderv(PH, "signalValue", -1)
PH <- PH[1:500]
ov <- list(PH18= PH$name,
           PH29= vl_intersectBed(PH, vl_importBed("db/peaks/cutnrun/PH_PH29_confident_peaks.narrowPeak"))$name,
           PHD9= vl_intersectBed(PH, vl_importBed("db/peaks/cutnrun/PH_PHD9_confident_peaks.narrowPeak"))$name,
           PHD11= vl_intersectBed(PH, vl_importBed("db/peaks/cutnrun/PH_PHD11_confident_peaks.narrowPeak"))$name)
vl_upset_plot(ov)

# K27 peaks
K27 <- vl_importBed("db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak")
setorderv(K27, "signalValue", -1)
K27 <- K27[1:100]
ov <- list(PH18= K27$name,
           PH29= vl_intersectBed(K27, vl_importBed("db/peaks/cutnrun/H3K27me3_PH29_confident_peaks.broadPeak")[signalValue>2])$name,
           PHD9= vl_intersectBed(K27, vl_importBed("db/peaks/cutnrun/H3K27me3_PHD11_confident_peaks.broadPeak")[signalValue>2])$name,
           PHD11= vl_intersectBed(K27, vl_importBed("db/peaks/cutnrun/H3K27me3_PHD9_confident_peaks.broadPeak")[signalValue>2])$name)
vl_upset_plot(ov)

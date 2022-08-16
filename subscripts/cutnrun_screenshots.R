setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

#-----------------------------------------------------------#
# Screenshots replicates
#-----------------------------------------------------------#
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
meta[, ChIP:= factor(ChIP, c("H3K27me3", "H2AK118Ub", "H3K27Ac", "H3K4me1", "H3K36me3", "IgG", "PH", "input"))]
meta[, cdition:= factor(cdition, c("PH18", "PHD11", "PHD9", "PH29"))]
meta[, col:= vl_palette_many_categ(.NGRP)[.GRP], ChIP]
setorderv(meta, c("ChIP", "cdition"))
# meta[ChIP=="H3K27me3", max:= 25]
# meta[ChIP=="H2AK118Ub", max:= 25]
# meta[ChIP=="H3K27Ac", max:= 20]
# meta[ChIP=="H3K4me1", max:= 20]
# meta[ChIP=="H3K36me3", max:= 20]
# meta[ChIP=="IgG", max:= 10]
# meta[ChIP=="PH", max:= 20]
# meta[ChIP=="input", max:= 10]

pdf("pdf/cutnrun/cutnrun_QC_screenshot_replicates.pdf", 
    width = 20, 
    height = 20)
par(mar= c(6, 10, 0, 2))
vl_screenshot(bed = GRanges("chr3R", IRanges(16.615e6, 17.025e6)),
              tracks = meta$bw, 
              col = meta$col,
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(6.625e6, 7.1e6)),
              meta$bw,
              col = meta$col,
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(5e6, 15e6)),
              meta$bw,
              col = meta$col,
              genome = "dm6")
dev.off()
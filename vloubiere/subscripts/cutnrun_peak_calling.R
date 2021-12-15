setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm6)

# Import data
dat <- data.table(file= list.files("db/bw/cutnrun_merge/", full.names = T))
dat[, cdition:= gsub("_merge.bw", "", basename(file))]
bins <- vl_binBSgenome(BSgenome.Dmelanogaster.UCSC.dm6, bin_size = 100) 

# Quantif bins and normalize
dat <- dat[, cbind(bins, vl_bw_coverage(bins, file)), .(cdition, file)]
dat <- na.omit(dat[score>0])
dat[, score:= log2(score)]
dat[, score:= score-median(score), cdition]

# Select enriched bins
K27Ac <- dat[score>1.5 & grepl("H3K27Ac", cdition)]
K27Ac <- vl_collapse_DT_ranges(K27Ac, mingap = 301)
K27Ac <- K27Ac[seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX") & end-start>190]
K27me3 <- dat[score>1.5 & grepl("H3K27me3", cdition)]
K27me3 <- vl_collapse_DT_ranges(K27me3, mingap = 10001)
K27me3 <- K27me3[seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX") & end-start>500]

pdf("pdf/cutnrun/screenshots_peak_calling_QC.pdf", width = 14)
vl_screenshot(GRanges("chr3R", IRanges(13.8e6, 14.2e6)),
              list.files("db/bw/cutnrun_merge/", "H3K27Ac", full.names= T),
              highlight_regions = K27Ac,
              genome = "dm6", 
              n_genes = 1)
vl_screenshot(GRanges("chr3R", IRanges(12.6e6, 17.04e6)),
              list.files("db/bw/cutnrun_merge/", "H3K27me3", full.names= T),
              highlight_regions = K27me3,
              genome = "dm6", 
              n_genes = 1)
dev.off()

export(K27Ac, "db/peaks/K27Ac_cutnrun.bed")
export(K27me3, "db/peaks/K27me3_cutnrun.bed")

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

bins <- vl_binBSgenome("dm6", bins_width = 500)
bins <- bins[seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")]
dat <- data.table(file= list.files("db/bw/cutnrun_reps_vl/", ".bw", full.names = T))
dat[, cdition:= gsub(".bw$","", basename(file))]
dat <- dat[, vl_bw_coverage(bins, file), .(cdition, file)]
dat[, bin:= rowid(cdition)]
mat <- dcast(dat, bin~cdition, value.var = "V1")
mat <- as.matrix(mat, 1) 
cor <- cor(mat)

pdf("pdf/cutnrun/PCC_new_cutnruns.pdf", 14, 13)
vl_heatmap(cor, display_numbers = T)
dev.off()

dir.create("db/peaks/PH_cutnrun", showWarnings = F)
PH18_rep1 <- "/mnt/f/_R_data/projects/epigenetic_cancer/db/sam/cutnrun/PH_PH18_rep1.sam"
PH18_rep2 <- "/mnt/f/_R_data/projects/epigenetic_cancer/db/sam/cutnrun/PH_PH18_rep2.sam"
PH29_rep1 <- "/mnt/f/_R_data/projects/epigenetic_cancer/db/sam/cutnrun/PH_PH29_rep1.sam"
PH29_rep2 <- "/mnt/f/_R_data/projects/epigenetic_cancer/db/sam/cutnrun/PH_PH29_rep2.sam"
system(paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 --outdir /mnt/d/_R_data/projects/epigenetic_cancer/db/peaks/PH_cutnrun/ -g dm --nomodel --keep-dup 1 --extsize 200 ",
              "-t ", PH18_rep1, " -c ", PH29_rep1, " -n ", " PH18_vs_PH29_rep1"))
system(paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 --outdir /mnt/d/_R_data/projects/epigenetic_cancer/db/peaks/PH_cutnrun/ -g dm --nomodel --keep-dup 1 --extsize 200 ",
              "-t ", PH18_rep2, " -c ", PH29_rep2, " -n ", " PH18_vs_PH29_rep2"))
system(paste0("/home/vloubiere/.local/bin/macs2 callpeak --keep-dup 1 --outdir /mnt/d/_R_data/projects/epigenetic_cancer/db/peaks/PH_cutnrun/ -g dm --nomodel --keep-dup 1 --extsize 200 ",
              "-t ", PH18_rep1, " ", PH18_rep2, " -c ", PH29_rep1, " ", PH29_rep2, " -n ", " PH18_vs_PH29_merge"))


peaks <- fread("db/peaks/PH_cutnrun/PH18_vs_PH29_merge_peaks.narrowPeak")
rep1 <- fread("db/peaks/PH_cutnrun/PH18_vs_PH29_rep1_peaks.narrowPeak")
rep2 <- fread("db/peaks/PH_cutnrun/PH18_vs_PH29_rep2_peaks.narrowPeak")
peaks$rep1 <- rep1[peaks, .N, .EACHI, on= c("V1", "V2<=V3", "V3>=V2")]$N>0
peaks$rep2 <- rep2[peaks, .N, .EACHI, on= c("V1", "V2<=V3", "V3>=V2")]$N>0
peaks <- peaks[(rep1 & rep2)]
fwrite(peaks[, V1:V10], 
       "db/peaks/PH_cutnrun/confident_PH_peaks_PH18_vs_PH29.bed", 
       sep= "\t",
       col.names = F, 
       row.names = F)

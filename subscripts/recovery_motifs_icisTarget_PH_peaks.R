# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]

#------------------------------#
# Compute and save regions
#------------------------------#
# PH
PH <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
PH[, start:= start+peak]
PH[, end:= start]
PH <- vl_closestBed(PH, dat)
PH <- PH[between(dist, -5000, 0) & !is.na(recovery.b), .(seqnames, start, end, strand, recovery= recovery.b)]
PH <- vl_resizeBed(PH, "center", 250, 250)
rtracklayer::export(GRanges(PH[recovery=="Recovery"]), 
                    "db/iCisTarget/Reversible_PH_peaks.bed")
rtracklayer::export(GRanges(PH[recovery=="noRecovery"]), 
                    "db/iCisTarget/Irreversible_PH_peaks.bed")

# Compare PH peaks
dat <- rbind(fread("db/iCisTarget/Irreversible_PH_peaks/statistics.tbl", sel= 1:8),
             fread("db/iCisTarget/Reversible_PH_peaks/statistics.tbl", sel= 1:8))
NES <- dcast(dat, 
             FeatureID+FeatureDescription+FeatureAnnotations~`#GeneSignatureID`, 
             value.var= "NES")
NES <- NES[FeatureAnnotations!="" & (Irreversible_PH_peaks>3 | Reversible_PH_peaks>3), 
           .SD[which.max(apply(.SD, 1, max))], FeatureAnnotations, 
           .SDcols= c("Irreversible_PH_peaks", "Reversible_PH_peaks")]

ggplot(NES, aes(Irreversible_PH_peaks, Reversible_PH_peaks, label = FeatureAnnotations)) + 
  geom_point(color = "red") + 
  geom_text_repel() +
  geom_abline(intercept = -1, slope = 1)+
  geom_abline(intercept = 1, slope = 1)+
  theme_classic()

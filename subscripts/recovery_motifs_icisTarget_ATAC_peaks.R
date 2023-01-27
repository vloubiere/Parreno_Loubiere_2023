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
# ATAC
RE <- vl_importBed("db/peaks/ATAC/ATAC_confident_peaks.narrowPeak")
RE[, start:= start+peak]
RE[, end:= start]
RE <- vl_closestBed(RE, dat)
RE <- RE[between(dist, -5000, 0) & !is.na(recovery.b), .(seqnames, start, end, strand, recovery= recovery.b)]
RE <- vl_resizeBed(RE, "center", 250, 250)
rtracklayer::export(GRanges(RE[recovery=="Recovery"]), 
                    "db/iCisTarget/Reversible_ATAC_peaks.bed")
rtracklayer::export(GRanges(RE[recovery=="noRecovery"]), 
                    "db/iCisTarget/Irreversible_ATAC_peaks.bed")

#------------------------------#
# Extract results
#------------------------------#
dat <- rbind(fread("db/iCisTarget/Irreversible_ATAC_peaks/statistics.tbl", sel= 1:8),
             fread("db/iCisTarget/Reversible_ATAC_peaks/statistics.tbl", sel= 1:8))
dat[grepl("^His2B:CG33878,His2B:CG33876,His2B:CG33874",FeatureAnnotations), FeatureAnnotations:= "His2B"]
dat[FeatureAnnotations=="Su(var)205,HP1e,HP1c,HP1b", FeatureAnnotations:= "Su(var)205/HP1b/c/e"]
dat[FeatureAnnotations=="Cfp1,CG17440,CG3347", FeatureAnnotations:= "Cfp1,CG3347,CG17440"]
dat[FeatureAnnotations=="Cnx14D,CG1924,Cnx99A", FeatureAnnotations:= "Cnx14D/99A,CG1924"]
dat[FeatureAnnotations=="Nf−YC,Nf−YB,Nf−YA", FeatureAnnotations:= "Nf−YA/B/C"]
dat[FeatureAnnotations=="Odc1,Odc2", FeatureAnnotations:= "Odc1/2"]
NES <- dcast(dat, 
             FeatureID+FeatureDescription+FeatureAnnotations~`#GeneSignatureID`, 
             value.var= "NES")
NES <- NES[FeatureAnnotations!="" & (Irreversible_ATAC_peaks>3 | Reversible_ATAC_peaks>3), 
           .SD[which.max(apply(.SD, 1, max))], FeatureAnnotations, 
           .SDcols= c("Irreversible_ATAC_peaks", "Reversible_ATAC_peaks")]

Cc <- fcase(NES$FeatureAnnotations %in% c("E(z)", "Su(z)12", "Adf1", "Trl"), "blue", 
            NES$FeatureAnnotations %in% c("zfh1", "Jra", "kay", "kay,Jra"), "red",
            default= "darkgrey")

pdf("pdf/motif_enrichment_iCisTarget_ATAC_peaks.pdf", 12, 12)
ggplot(NES, aes(Irreversible_ATAC_peaks, Reversible_ATAC_peaks, label = FeatureAnnotations)) + 
  xlim(-1,7) +
  ylim(-1,7) +
  geom_point(color = Cc, size= 2) + 
  geom_text_repel(max.overlaps= 20, col= Cc, size= ifelse(Cc=="darkgrey", 6, 7.5)) +
  geom_abline(intercept = -1, slope = 1, linetype="dashed")+
  geom_abline(intercept = 1, slope = 1, linetype="dashed")+
  theme_classic()
dev.off()
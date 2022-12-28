# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")[, .(seqnames, start, end, strand, recovery, FPKM_PH18)]
# Select recovery genes whose TSS is bound by PRC1
dat <- vl_resizeBed(dat[!is.na(recovery)], "start", 250, 250)
PRC1 <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
dat <- vl_intersectBed(dat, PRC1)
dat[, recovery:= switch(recovery, 
                        "noRecovery"= "Irreversible",
                        "Recovery"= "Reversible"), recovery]

Cc <- c("darkorchid2", "chartreuse3")

pdf("pdf/Figure_4_boxplots_FPKM_PHbound_genes.pdf", 
    width = 2.5, 
    height = 3)
par(mar= c(4,5,2,2),
    las= 1,
    tcl= -0.2,
    mgp= c(1.5, 0.5, 0))
vl_boxplot(FPKM_PH18~recovery, 
           dat,
           compute_pval= list(c(1,2)),
           col= adjustcolor(Cc, 0.6),
           tilt.names= T,
           srt= 30,
           ylab= "RNA-Seq FPKM")
dev.off()
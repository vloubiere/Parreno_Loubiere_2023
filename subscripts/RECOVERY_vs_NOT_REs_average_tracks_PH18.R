setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

#------------------------------------------------#
# Import data
#------------------------------------------------#
# Import and resize assigned ATAC-Seq peaks
ATAC <- readRDS("Rdata/motif_counts_ATAC_sites.rds")[group=="ATAC", c(1,3), with=F]
ATAC[, c("seqnames", "start", "end"):= as.data.table(GRanges(coor))[, .(seqnames, start, end)]]
ATAC[, start:= round(rowMeans(.SD)-1000), .SDcols= c("start", "end")]
ATAC[, end:= start+2000]
ATAC$coor <- NULL
# make dat object
dat <- data.table(track= c("db/bw/SA_2020/PC_ED_merge.bw",
                           "db/bw/SA_2020/PSC_ED_merge.bw",
                           "db/bw/SA_2020/PH_ED_merge.bw",
                           "db/bw/cutnrun/PH_PH18_rep1.bw",
                           "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
                           "db/bw/cutnrun/H2AK118Ub_PH18_merge.bw",
                           "db/bw/cutnrun/H3K27me3_PH18_merge.bw",
                           "db/bw/cutnrun/H3K36me3_PH18_merge.bw",
                           "db/bw/cutnrun/H3K4me1_PH18_merge.bw",
                           "db/bw/cutnrun_EcR/EcR_-6hAPF_merge.bw",
                           "db/bw/cutnrun_EcR/EcR_+6hAPF_merge.bw"))
dat <- dat[, fread("Rdata/RECOVERY_NORECOVERY_genes.txt"), (dat)]
dat <- dat[, .(track, FBgn, PHD9_RECOVERY, PHD11_RECOVERY)]
dat <- ATAC[dat, on= "FBgn"]
dat <- melt(dat, 
            measure.vars = c("PHD9_RECOVERY", "PHD11_RECOVERY"), 
            value.name = "class")
dat[, cdition:= gsub("_merge.bw|_rep1.bw", "", basename(track)), track]
# Quantif tracks
dat[, quantif:= {
  vl_bw_coverage(data.table(seqnames, start, end), track)
}, track]
dat <-na.omit(dat)

#-----------------------------------------#
# Plot
#-----------------------------------------#
pdf("pdf/Figures/PH18_CUTNRUN_enrich_REs_revert_vs_not.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
dat[, {
  # Plot average track for each class
  .q <- vl_bw_average_track(bed= data.table(seqnames, start, end), 
                            tracks= track,
                            set_IDs = class,
                            upstream = 2500,
                            downstream = 2500,
                            stranded = T,
                            center_label = "Center", 
                            legend.cex = 0.6, 
                            legend= F)
  # Legend
  leg <- .q[, .(N= length(unique(region_ID))), .(set_IDs, col)]
  leg[, legend("topright",
               legend = paste0(ifelse(set_IDs, "RECOVERY (", "NO RECOVERY ("), N, ")"),
               fill= col,
               bty= "n")]
  title(main= paste(cdition, variable))
  # Boxplot
  x <- split(quantif, class)
  vl_boxplot(x,
             boxcol= vl_palette_few_categ(2),
             ylab= "Enrichment", 
             compute_pval= list(c(1,2)))
  print(".")
}, .(variable, track, cdition)]
dev.off()

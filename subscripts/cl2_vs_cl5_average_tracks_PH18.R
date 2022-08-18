setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Make dat object
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery), .(recovery, seqnames, start, end, strand)]
dat <- dat[, .(track= c("db/bw/SA_2020/PC_ED_merge.bw",
                        "db/bw/SA_2020/PSC_ED_merge.bw",
                        "db/bw/SA_2020/PH_ED_merge.bw",
                        "db/bw/cutnrun/PH_PH18_rep1.bw",
                        "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
                        "db/bw/cutnrun/H2AK118Ub_PH18_merge.bw",
                        "db/bw/cutnrun/H3K27me3_PH18_merge.bw",
                        "db/bw/cutnrun/H3K36me3_PH18_merge.bw",
                        "db/bw/cutnrun/H3K4me1_PH18_merge.bw",
                        "db/bw/cutnrun_EcR/EcR_-6hAPF_merge.bw",
                        "db/bw/cutnrun_EcR/EcR_+6hAPF_merge.bw")), (dat)]
dat[, cdition:= gsub("_merge.bw|_rep1.bw", "", basename(track)), track]
dat[, protein:= grepl("^PC|^PH|^PSC|^EcR", cdition)]
# Quantif tracks
dat[, quantif:= {
  if(protein)
  {
    upstream <- 250
    downstream <- 250
  }else
  {
    upstream <- 2500
    downstream <- 10000
  }
  vl_bw_coverage(vl_resizeBed(data.table(seqnames, start, end, strand), 
                              "start", 
                              upstream = upstream, 
                              downstream = downstream), 
                 track)
}, .(track, protein)]

#-----------------------------------------#
# Plot
#-----------------------------------------#
pdf("pdf/Figures/PH18_CUTNRUN_enrich_cl2_vs_cl5.pdf",
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
  .q <- vl_bw_average_track(bed= data.table(seqnames, start, end= start, strand), 
                            tracks= track,
                            set_IDs = recovery,
                            upstream = 2500,
                            downstream = 10000,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            legend= F)
  # Legend
  leg <- .q[, .(N= length(unique(region_ID))), .(set_IDs, col)]
  leg[, legend("topright",
               legend = paste0(ifelse(set_IDs, "RECOVERY (", "NO RECOVERY ("), N, ")"),
               fill= col,
               bty= "n")]
  title(main= cdition)
  # Boxplot
  x <- split(quantif, recovery)
  vl_boxplot(x,
             boxcol= vl_palette_few_categ(2),
             ylab= "Enrichment", 
             compute_pval= list(c(1,2)),
             notch= T)
  print(".")
}, .(track, cdition)]
dev.off()

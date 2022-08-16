setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Make dat object
dat <- data.table(track= c("db/bw/cutnrun/H2AK118Ub_PH18_merge.bw",
                           "db/bw/cutnrun/H2AK118Ub_PHD11_merge.bw",
                           "db/bw/cutnrun/H2AK118Ub_PH29_merge.bw",
                           "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
                           "db/bw/cutnrun/H3K27Ac_PHD11_merge.bw",
                           "db/bw/cutnrun/H3K27Ac_PH29_merge.bw",
                           "db/bw/cutnrun/H3K27me3_PH18_merge.bw",
                           "db/bw/cutnrun/H3K27me3_PHD11_merge.bw",
                           "db/bw/cutnrun/H3K27me3_PH29_merge.bw",
                           "db/bw/cutnrun/H3K36me3_PH18_merge.bw",
                           "db/bw/cutnrun/H3K36me3_PHD11_merge.bw",
                           "db/bw/cutnrun/H3K36me3_PH29_merge.bw",
                           "db/bw/cutnrun/H3K4me1_PH18_merge.bw",
                           "db/bw/cutnrun/H3K4me1_PHD11_merge.bw",
                           "db/bw/cutnrun/H3K4me1_PH29_merge.bw"),
                  ChIP= rep(c("H2AK118Ub", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1"), each= 3),
                  cdition= rep(c("PH18", "PHD11", "PH29"), 5),
                  ymin= rep(c(1.5,1,2,0.5,1.5), each= 3),
                  ymax= rep(c(7,6,20,4.5,10.5), each= 3))
dat <- dat[, fread("Rdata/cl2_cl5_RECOVERY_genes.txt"), (dat)]
# Add prom coor
genes <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[type=="gene", .(FBgn= gene_id, seqnames, start, end, strand, symbol= gene_symbol)]
proms <- vl_resizeBed(genes, "start", 0, 0)
dat <- proms[dat, on= "FBgn"]
# Quantif tracks
dat[, quantif:= {
  coor <- vl_resizeBed(data.table(seqnames, start, end, strand), 
                       "center", 
                       upstream = 2500, 
                       downstream = 10000)
  vl_bw_coverage(coor, track)
}, track]

#-----------------------------------------#
# Plot
#-----------------------------------------#
pdf("pdf/Figures/compare_CUTNRUN_cditions_cl2_vs_cl5.pdf",
    height = 2.75, 
    width = 9)
layout(matrix(1:4, ncol=4, byrow = T), 
       widths= rep(c(1,0.5), 2))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
dat[, {
  # Plot average track for each class
  .q <- vl_bw_average_track(bed= data.table(seqnames, start, end), 
                            tracks= unique(track),
                            upstream = 2500,
                            downstream = 10000,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            legend= F,
                            ylim= c(ymin, ymax))
  # Legend
  .q[, cdition:= tstrsplit(basename(file), "_", keep = 2)]
  leg <- unique(.q[, .(cdition, col)])
  leg[, legend("topright",
               legend = cdition,
               fill= col,
               bty= "n")]
  title(main= paste(ChIP, 
                    ifelse(RECOVERY, "RECOVERY", "NO RECOVERY")))
  # Boxplot
  x <- split(quantif, cdition)
  vl_boxplot(x,
             boxcol= vl_palette_few_categ(3), 
             tilt.names= T,
             ylab= "Enrichment", 
             compute_pval= list(c(1,2), c(1,3)))
  print(".")
}, keyby= .(ChIP, RECOVERY, ymin, ymax)]
dev.off()

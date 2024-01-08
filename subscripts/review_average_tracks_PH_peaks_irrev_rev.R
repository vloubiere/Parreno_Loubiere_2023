setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import K27 domains ----
K27 <- readRDS("Rdata/K27_domains_classif_reversible_irreversible.rds")

# Import PH peaks that overlap K27 domains ----
PH <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
PH[, start:= start+peak]
PH[, end:= start]
PH <- vl_intersectBed(PH, K27[reversible+irreversible>0])

# Import TSSs ----
gene <- readRDS("Rdata/final_gene_features_table.rds")
gene <- gene[class %in% c("Irreversible", "Reversible")]
TSS <- vl_resizeBed(gene, "start", 0, 0)

# Assign to closest TSS ----
peaks <- vl_closestBed(PH, TSS)
# Make sure that each peak has a unique assignment
if(nrow(unique(peaks[, .(seqnames, start, end, class.b)])) != nrow(peaks[, .(seqnames, start, end)]))
  stop("Some PH peaks have dual assignment")
peaks <- peaks[abs(dist)<25000]

# Quantification window ----
peaks <- vl_resizeBed(peaks, "center", 250, 250)
peaks[, class:= paste0(class.b, " (", .N, ")"), class.b]
peaks[, class:= factor(class, sort(unique(class)))]

# Import tracks ----
tracks <- fread("db/bw/list_files.txt")[ChIP=="PH"]

# Plots ----
Cc <- c("rosybrown1", "palegreen3")
ncol <- 1
nrow <- 3

pdf("pdf/review_average_tracks_PH_peaks.pdf",
    2.35*ncol,
    1.75*nrow)
layout(matrix(1:(nrow*ncol*2), ncol= ncol*2, byrow = T),
       widths = rep(c(1, .6), ncol))
par(las= 1,
    mgp= c(1, .3, 0),
    tcl= -0.1,
    mai= c(.5, .4, .5, .1),
    cex.axis= 7/12,
    cex.lab= 9/12)
tracks[, {
  # browser()
  vl_bw_average_track(peaks,
                      set.IDs = peaks$class,
                      tracks = bw_merge,
                      upstream = 2000,
                      downstream = 2000,
                      legend.pos = "topright",
                      legend.cex = .6,
                      center.label = "PH peaks\nsummit",
                      ylab= "PH average signal",
                      col= Cc)
  abline(v= c(-250,250), lty= "11")
  title(main= title)
  .c <- vl_bw_coverage(peaks, bw_merge)
  .c <- split(.c, peaks$class)
  vl_boxplot(.c,
             col= adjustcolor(Cc, .6),
             tilt.names = T, 
             compute.pval= list(c(1,2)),
             ylab= "Enrichment",
             lwd= 0.75)
  print(title)
}, .(bw_merge, title)]
dev.off()
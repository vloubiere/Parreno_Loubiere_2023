setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import genes ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[!is.na(class)]
setorderv(dat, "class")
dat[, class:= paste0(class, " (", formatC(.N, big.mark= ","), ")"), class]
dat[, class:= factor(class, unique(class))]

# Quantif ranges ----
TSS <- vl_resizeBed(dat, "start", 0, 0)
promExt <- c(750, 250)
prom <- vl_resizeBed(dat, "start", promExt[1], promExt[2])
htmExt <- c(2500, 2500)
ext <- vl_resizeBed(dat, "start", htmExt[1], htmExt[2])

# Import tracks ----
tracks <- fread("db/bw/list_files.txt")
tracks <- tracks[ChIP %in% c("PH", "H3K27me3", "H2AK118Ub", "H3K27Ac")]
tracks[, ChIP:= factor(ChIP, c("PH", "H3K27me3", "H2AK118Ub", "H3K27Ac"))]
setorderv(tracks, "ChIP")
tracks[, broad:= ChIP!="PH"]

# Plots ----
Cc <- c("lightgrey", "rosybrown1", "palegreen3")
ncol <- 3
nrow <- 4

pdf("pdf/review_average_tracks_promoters.pdf",
    2.5*ncol,
    1.35*nrow)
layout(matrix(1:(nrow*ncol*2), ncol= ncol*2, byrow = T),
       widths = rep(c(1, .35), ncol))
par(las= 1,
    mgp= c(1.2, .3, 0),
    tcl= -0.1,
    cex.axis= 7/12,
    cex.lab= 9/12)
tracks[, {
  # Average track ----
  par(mai= c(.425, .275, .2, .015))
  vl_bw_average_track(TSS,
                      set.IDs = dat$class,
                      tracks = bw_merge,
                      upstream = ifelse(broad, 10000, 2500),
                      downstream = ifelse(broad, 25000, 10000),
                      legend = .GRP==1,
                      legend.pos = "topright",
                      legend.cex = .6,
                      ignore.strand = F,
                      center.label = "TSS",
                      ylab= "Average signal",
                      col= Cc)
  title(main= title,
        line= .5)
  if(broad)
  {
    .c <- vl_bw_coverage(ext, bw_merge)
    abline(v= htmExt*c(-1,1), lty= "11")
  }else
  {
    .c <- vl_bw_coverage(prom, bw_merge)
    abline(v= promExt*c(-1,1), lty= "11")
  }
  
  # Boxplot quantif ----
  par(mai= c(.425, .275, .2, .035))
  .c <- split(.c, dat$class)
  vl_boxplot(.c,
             col= adjustcolor(Cc, .6),
             tilt.names = T, 
             compute.pval= list(c(2, 3)),
             ylab= "Enrichment",
             lwd= 0.75)
  print(title)
}, .(bw_merge, title, broad)]
dev.off()
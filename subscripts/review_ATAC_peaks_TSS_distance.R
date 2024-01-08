setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import peaks and TSSs ----
peaks <- readRDS("Rdata/final_ATAC_table.rds")
genes <- readRDS("Rdata/final_gene_features_table.rds")
TSS <- vl_resizeBed(genes, "start", 0, 0)

# Assign peaks to closest genes ----
dat <- vl_closestBed(peaks, TSS)
res <- dat[, sum(abs(dist)>1000, na.rm= T)/.N*100, keyby= cl]

# Cc <- c("darkgrey", "khaki3", "rosybrown1", "palegreen3", "darkgreen", "cornflowerblue", "black")
Cc <- c("darkgrey", "rosybrown1", "palegreen3", "cornflowerblue")

# Plot ----
pdf("pdf/review_ATAC_peaks_distance_to_TSS.pdf", 3, 3)
vl_par(mai= c(.9, 1.1, .9, 1.1))
bar <- barplot(res$V1,
        col= adjustcolor(Cc, .6),
        ylab= "Percentage of distal REs (>1kb from closest TSS)",
        border= NA)
vl_tilt_xaxis(bar,
              labels= res$cl)
dev.off()
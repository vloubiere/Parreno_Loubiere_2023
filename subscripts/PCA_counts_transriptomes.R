setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(DESeq2)
require(data.table)
require(vlfunctions)

# Import data ----
dat <- readRDS("db/dds/RNA/epiCancer_noGFP_dds.rds")
dat <- as.data.table(counts(dat, normalized= T))
mat <- t(log2(dat+1))
# PCA
pca <- as.data.table(prcomp(mat)$x, keep.rownames = "cdition")
cditions <- gsub("_rep1|_rep2|_rep3", "", pca$cdition)
cditions <- factor(cditions, levels= unique(cditions))
Cc <- rainbow(length(levels(cditions)))[as.numeric(cditions)]

# Plot ----
pdf("pdf/RNA_PCA_counts.pdf",
    width = 6.5,
    height = 5.5)
par(las= 1,
    mai= c(1.02, 0.82, 0.82, 1.92))
plot(pca$PC1, 
     pca$PC2, 
     xlab= "PC1",
     ylab= "PC2",
     col= adjustcolor(Cc, 0.6), 
     pch= 16,
     cex= 2)
legend(par("usr")[2],
       par("usr")[4],
       levels(cditions),
       col= adjustcolor(unique(Cc), 0.6),
       pch= 16,
       bty= "n",
       pt.cex= 2,
       xpd= T)
dev.off()
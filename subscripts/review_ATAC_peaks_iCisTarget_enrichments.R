setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import ATAC peaks ----
ATAC <- readRDS("Rdata/final_ATAC_table.rds")
ATAC <- vl_resizeBed(ATAC,
                     "center",
                     upstream = 250,
                     downstream = 250,
                     genome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)

# Export peaks ----
ATAC[, file:= paste0("db/iCisTarget/review/ATAC_peaks/", gsub(" ", "_", cl), "_peaks.bed"), cl]
ATAC[, {
  fwrite(.SD[, .(seqnames, start, end)],
         file,
         col.names = F,
         row.names = F,
         sep= "\t",
         quote= F)
}, file]

# Download results ----
# vl_download_iCistarget(url = "https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/reports/89c22109d9c7989efd3c6dfb16be4c09c79aedfe/archive.zip",
#                        dir = "db/iCisTarget/review/ATAC_peaks/")
# vl_download_iCistarget(url = "https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/reports/c14aaa0c98758571713b6c80a3ac4461a229dbdc/archive.zip",
#                        dir = "db/iCisTarget/review/ATAC_peaks/")
# vl_download_iCistarget(url = "https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/reports/45517214f7645e8809a663c4371248b5d8312b25/archive.zip",
#                        dir = "db/iCisTarget/review/ATAC_peaks/")
# vl_download_iCistarget(url = "https://gbiomed.kuleuven.be/apps/lcb/i-cisTarget/reports/8b1658ce4d635041d74c41cbbdc0599325db9c52/archive.zip",
#                        dir = "db/iCisTarget/review/ATAC_peaks/")


# make final table ----
if(!file.exists("db/iCisTarget/review/ATAC_peaks/results.rds"))
{
  res <- data.table(file= list.files("db/iCisTarget/review/ATAC_peaks/", "statistics.tbl", recursive= T, full.names = T))
  res[, variable:= tstrsplit(file, "/", keep= 6)]
  res[, variable:= gsub("_genes", "", variable)]
  res[, variable:= gsub("^cluster", "Cluster", variable)]
  res <- res[, fread(file, sel= 1:8), .(variable, file)]
  setnames(res, "#GeneSignatureID", "cdition")
  res[, cdition:= gsub("^cluster_|_ATAC_peaks$|_peaks$", "", cdition)]
  res <- res[FeatureDatabase=="PWMs"]
  # Get names and simplify
  res[, name:= ifelse(FeatureAnnotations=="", FeatureDescription, FeatureAnnotations)]
  res[grepl("[", name, fixed = T), name:= tstrsplit(name, split= "\\[", perl = T, keep= 1)]
  res[, name:= gsub(" ", "", name)]
  # Expressed genes
  symb <- readRDS("Rdata/final_gene_features_table.rds")
  symb <- symb[, symbol[apply(.SD, 1, function(x) any(x>1))], .SDcols= patterns("^FPKM")] #Symbol expressed genes
  res <- res[name %in% symb]
  saveRDS(res,
          "db/iCisTarget/review/ATAC_peaks/results.rds")
}

# Import data and select motifs for heatmap ----
dat <- readRDS("db/iCisTarget/review/ATAC_peaks/results.rds")
dat[, variable:= factor(variable, c("Irreversible", "Reversible", "Decreased"))]
dat <- na.omit(dat)
# Select top motifs
dat <- dat[FeatureID %in% dat[Rank<50 & NES>5.5, FeatureID]]
mat <- dcast(dat,
             name~variable,
             value.var = "NES",
             fun.aggregate = max)
mat <- as.matrix(mat, 1)
mat[is.na(mat)] <- 0
num <- round(mat, 1)
num[num<3] <- ""

# TFs transcriptional levels ----
FPKM <- readRDS("Rdata/final_gene_features_table.rds")
FPKM <- FPKM[, .(symbol, FPKM_PH18, FPKM_PH29, FPKM_PHD11)]

# Import data and select motifs for scatterplot ----
dat <- readRDS("db/iCisTarget/review/ATAC_peaks/results.rds")
dat <- dat[variable %in% c("Irreversible", "Reversible")]
dat <- dat[FeatureID %in% dat[NES>4, FeatureID]]
pl <- dcast(dat,
             name~variable,
             value.var = "NES",
            fun.aggregate = max,
            fill = 0)

# Plot ----
pdf("pdf/review_ATAC_peaks_iCisTarget_enrichments.pdf", 6, 6)
par(las= 2,
    mgp= c(1.5, 0.35, 0),
    mai= c(2,3.2,.5,2),
    cex.axis= 9/12)
hm <- vl_heatmap(mat,
                 cluster.rows = T,
                 cluster.cols= F,
                 col= c("white", "white", "tomato"),
                 breaks= c(0, 2, 7),
                 display.numbers = T,
                 display.numbers.matrix = num,
                 display.numbers.cex = .6,
                 legend.cex = 6/12,
                 legend.title = "NES")
act <- FPKM[hm$rows[(order)], symbol:FPKM_PHD11, on= "symbol==name"]
logAct <- log2(as.matrix(act, 1)+1)
# Heatmap FPKM
vl_heatmap(logAct,
           cluster.rows= F,
           cluster.cols= F,
           col= c("blue", "yellow"),
           breaks = c(0,7),
           display.numbers= T,
           display.numbers.matrix = round(logAct, 1),
           display.numbers.cex = .6)
par(mai= rep(2, 4),
    las= 1,
    tcl= -0.1)
pl[, {
  vl_repelScatterplot(Irreversible,
                      Reversible,
                      labels = name,
                      label.cex = 5/12,
                      pch= 16,
                      col= adjustcolor(ifelse(name %in% c("Stat92E", "zfh1", "Abd-B", "cad"), "tomato", "lightgrey"), .8),
                      xlim= c(-0.5, 8),
                      ylim= c(-0.5, 8),
                      point_size = .1)
  abline(1, 1, lty= "11")
  abline(-1, 1, lty= "11")
}]
dev.off()
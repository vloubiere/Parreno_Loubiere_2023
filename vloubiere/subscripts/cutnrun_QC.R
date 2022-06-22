setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)
require(kohonen)

# SPIKE IN fraction
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
meta <- meta[ChIP %in% c("H3K27me3", "H3K27Ac", "IgG") & Comment!="failed" & Suffix==""]
meta[, ChIP_reads:= fread(cmd= paste0("wc -l ", ChIP_bed))$V1, ChIP_bed]
meta[, spikein_reads:= fread(cmd= paste0("wc -l ", spikein_bed))$V1, spikein_bed]
pl <- unique(meta[, .(name= paste0(ChIP, "_", cdition, "_", rep), ChIP_reads, spikein_reads)])
pl <- t(as.matrix(pl, 1))

pdf("pdf/cutnrun/spikein_percentage.pdf")
layout(matrix(c(1,3,2), nrow= 3), 
       heights = c(1.7,0.2,2))
Cc <- c("grey20", "grey80")
par(mar= c(2,5,2,2))
barplot(pl, 
        beside = T,
        las= 2, 
        col= Cc,
        ylab= "N reads\n", 
        xaxt= "n")
legend("topright",
       fill= Cc,
       legend= c("ChIP", "Spike-in"),
       bty= "n",
       xpd= T)
par(mar= c(12,5,2,2))
bar <- barplot(pl["spikein_reads",]/pl["ChIP_reads",]*100, 
               ylim= c(0, 2.5),
               beside = T,
               las= 2,
               ylab= "% spike-in reads")
text(bar, 
     grconvertY(1, "nfc", "user"),
     paste0(round(pl["spikein_reads",]/pl["ChIP_reads",]*100, 2), "%"),
     pos= 1,
     xpd= T)
dev.off()

#-----------------------------------------------------------#
# Correlations
#-----------------------------------------------------------#
bins <- vl_binBSgenome(BSgenome.Dmelanogaster.UCSC.dm6, bins_width = 5000) 
dat <- meta[, SJ(idx= seq(nrow(bins)), 
                 score= vl_bw_coverage(bins, bw)), .(bw, cdition= paste0(ChIP, "_", cdition, "_", rep))]
dat <- dcast(dat, 
             idx~cdition, 
             value.var = "score")
cols <- grep("rep", names(dat), value = T)

# Heatmap
pdf("pdf/cutnrun/PCC_cutnrun_replicates.pdf", 
    width= 9,
    height= 9)
vl_heatmap(cor(na.omit(dat[, ..cols])), 
           display_numbers = T, 
           col= c("blue", "yellow"),
           legend_title = "PCC")
dev.off()

#-----------------------------------------------------------#
# PCA
#-----------------------------------------------------------#
cols <- cols[!grepl("^IgG", cols)]
pca <- as.data.table(prcomp(na.omit(dat[, ..cols]))$rotation,
                     keep.rownames= T)

pca[grep("18", rn), Cc:= "grey50"]
pca[grep("D9", rn), Cc:= "gold"]
pca[grep("D11", rn), Cc:= "limegreen"]
pca[grep("29", rn), Cc:= "orangered"]
pca[grep("ChIP-Seq", rn), Cc:= "white"]
pca[grep("H3K27me3", rn), Cc:= colorRampPalette(c(Cc, "blue"))(3)[2], Cc]
pca[grep("H3K27Ac", rn), Cc:= colorRampPalette(c(Cc, "red"))(3)[2], Cc]

pdf("pdf/cutnrun/PCA_cutnrun_replicates.pdf", 
    width = 11)
par(mar= c(4,4,1,20))
plot(pca[, PC1],
     pca[, PC2], 
     col= pca[, Cc],
     pch= 16,
     cex= 2,
     xlab= "PC1",
     ylab= "PC2")
legend(par("usr")[2],
       par("usr")[4],
       xpd= T,
       pch= 16,
       cex= 1.5,
       legend = unique(gsub("_rep1|_rep2", "", pca[, rn])),
       col= unique(pca[, Cc]),
       bty= "n")
dev.off()


#-----------------------------------------------------------#
# Screenshots replicates
#-----------------------------------------------------------#
meta[, ChIP:= factor(ChIP, c("H3K27Ac", "H3K27me3", "IgG"))]
meta[, cdition:= factor(cdition, c("PH18", "PHD11", "PHD9", "PH29"))]
setorderv(meta, c("ChIP", "cdition"))
meta[ChIP=="H3K27Ac", max:= 20]
meta[ChIP=="H3K27me3", max:= 25]
meta[ChIP=="IgG", max:= 10]

pdf("pdf/cutnrun/cutnrun_QC_screenshot_replicates.pdf", 
    width = 10, 
    height = 10)
par(mar= c(6, 10, 0, 2))
vl_screenshot(bed = GRanges("chr3R", IRanges(16.615e6, 17.025e6)),
              tracks = meta$bw, 
              max= meta$max,
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(6.625e6, 7.1e6)),
              meta$bw,
              max= meta$max,
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(5e6, 15e6)),
              meta$bw,
              max= meta$max,
              genome = "dm6")
dev.off()

#-----------------------------------------------------------#
# Screenshots merge
#-----------------------------------------------------------#
tracks <- unique(meta[, .(bw_merge, max)])

pdf("pdf/cutnrun/cutnrun_QC_screenshot_merged_tracks.pdf", 
    width = 10, 
    height = 6)
par(mar= c(6, 10, 0, 1))
vl_screenshot(GRanges("chr3R", IRanges(16.615e6, 17.025e6)),
              tracks$bw_merge, 
              max= tracks$max,
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(6.625e6, 7.1e6)),
              tracks$bw_merge,
              max= tracks$max,
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(5e6, 15e6)),
              tracks$bw_merge,
              max= tracks$max,
              genome = "dm6")
dev.off()

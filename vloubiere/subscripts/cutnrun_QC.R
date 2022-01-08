setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)
require(kohonen)

# SPIKE IN fraction
meta <- fread("Rdata/metadata_cutnrun_final.txt")
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
bar <- barplot(pl, 
               ylim= c(0, 5e4),
               beside = T,
               las= 2, 
               col= Cc,
               ylab= "N reads\n")
rect(grconvertX(0, "nfc", "user"),
     grconvertY(0.95, "nfc", "user"),
     grconvertX(1, "nfc", "user"),
     grconvertY(1, "nfc", "user"), 
     col = "white", 
     border= NA,
     xpd= T)
text(colMeans(bar), 
     grconvertY(1, "nfc", "user"),
     paste0(round(pl["spikein_reads",]/pl["ChIP_reads",]*100, 2), "%"),
     pos= 1,
     xpd= T)
dev.off()

#-----------------------------------------------------------#
# Correlations
#-----------------------------------------------------------#
# Import data with replicates
if(!file.exists("Rdata/gw_250bp_bins_quantif_cutnrun_replicates.rds"))
{
  dat <- data.table(file= c("db/bw/SA_2020/H3K27Ac_ED_merge.bw",
                            "db/bw/SA_2020/H3K27me3_ED_merge.bw",
                            list.files("db/bw/cutnrun_reps_vl/", full.names = T)))
  dat[, cdition:= gsub(".bw", "", basename(file))]
  dat[, cdition:= gsub("merge", "ChIP-Seq", cdition)]
  bins <- vl_binBSgenome(BSgenome.Dmelanogaster.UCSC.dm6, bin_size = 250) 
  dat <- dat[, cbind(bins, score= vl_bw_coverage(bins, file)), .(cdition, file)]
  dat <- dcast(dat, 
               seqnames+start+end~cdition, 
               value.var = "score")
  saveRDS(dat, 
          "Rdata/gw_250bp_bins_quantif_cutnrun_replicates.rds")
}else if(!exists("dat"))
  dat <- readRDS("Rdata/gw_250bp_bins_quantif_cutnrun_replicates.rds")
dat <- na.omit(dat)

# Heatmap
cols <- grep("rep", names(dat), value = T)

pdf("pdf/cutnrun/PCC_cutnrun_replicates_wo_ChIP-Seq.pdf")
vl_heatmap(cor(na.omit(dat[, ..cols])), 
           display_numbers = T, 
           col= c("blue", "yellow"),
           legend_title = "PCC")
dev.off()

pdf("pdf/cutnrun/PCC_cutnrun_replicates_w_ChIP-Seq.pdf")
vl_heatmap(cor(na.omit(dat[, `H3K27Ac_ED_ChIP-Seq`:H3K27me3_PHD9_rep2])), 
           display_numbers = T, 
           col= c("blue", "yellow"),
           legend_title = "PCC")
dev.off()

#-----------------------------------------------------------#
# PCA
#-----------------------------------------------------------#
# Without ChIP-Seq
pca <- as.data.table(prcomp(dat[, ..cols])$rotation,
                     keep.rownames= T)
pca[grep(18, rn), Cc:= "grey50"]
pca[grep("D9", rn), Cc:= "gold"]
pca[grep("D11", rn), Cc:= "limegreen"]
pca[grep("29", rn), Cc:= "orangered"]
pca[grep("ChIP-Seq", rn), Cc:= "white"]
pca[grep("H3K27me3", rn), Cc:= colorRampPalette(c(Cc, "blue"))(3)[2], Cc]
pca[grep("H3K27Ac", rn), Cc:= colorRampPalette(c(Cc, "red"))(3)[2], Cc]

pdf("pdf/cutnrun/PCA_cutnrun_replicates_wo_ChIP-Seq.pdf", 
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
       legend = pca[, rn],
       col= pca[, Cc],
       bty= "n")
dev.off()


#-----------------------------------------------------------#
# Screenshots replicates
#-----------------------------------------------------------#
tracks <- data.table(file= c("db/bw/SA_2020/PH_ED_merge.bw",
                             "db/bw/SA_2020/EZ_CNSID_BL_merge.bw",
                             "db/bw/SA_2020/SUZ12_ED_merge.bw",
                             list.files("db/bw/cutnrun_reps_vl/", full.names = T)))
pdf("pdf/cutnrun/screenshot_replicates.pdf", 
    width = 10, 
    height = 10)
vl_screenshot(GRanges("chr3R", IRanges(16.615e6, 17.025e6)),
              tracks$file, 
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(6.625e6, 7.1e6)),
              tracks$file, 
              genome = "dm6")
vl_screenshot(GRanges("chr3R", IRanges(5e6, 15e6)),
              tracks$file, 
              genome = "dm6")
dev.off()

#-----------------------------------------------------------#
# Screenshots merge
#-----------------------------------------------------------#
tracks <- data.table(file= c("db/bw/SA_2020/PH_ED_merge.bw",
                             "db/bw/SA_2020/EZ_CNSID_BL_merge.bw",
                             "db/bw/SA_2020/SUZ12_ED_merge.bw",
                             list.files("db/bw/cutnrun_merge_vl/", full.names = T)))
pdf("pdf/cutnrun/screenshot_merge.pdf", 
    width = 10, 
    height = 6)
vl_screenshot(GRanges("chr3R", IRanges(16.615e6, 17.025e6)),
              tracks$file, 
              genome = "dm6", 
              n_genes = 1)
vl_screenshot(GRanges("chr3R", IRanges(6.625e6, 7.1e6)),
              tracks$file, 
              genome = "dm6",
              n_genes = 1)
vl_screenshot(GRanges("chr3R", IRanges(5e6, 15e6)),
              tracks$file, 
              genome = "dm6",
              n_genes = 1)
dev.off()

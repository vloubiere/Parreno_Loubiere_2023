setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)
require(kohonen)

# Import data with replicates
if(!file.exists("Rdata/gw_250bp_bins_quantif_cutnrun_replicates.rds"))
{
  dat <- data.table(file= c("db/bw/SA_2020/H3K27Ac_ED_merge.bw",
                            "db/bw/SA_2020/H3K27me3_ED_merge.bw",
                            list.files("db/bw/cutnrun_reps/", full.names = T)))
  dat[, cdition:= gsub(".bw", "", basename(file))]
  dat[, cdition:= gsub("merge", "ChIP-Seq", cdition)]
  bins <- vl_binBSgenome(BSgenome.Dmelanogaster.UCSC.dm6, bin_size = 250) 
  dat <- dat[, cbind(bins, vl_bw_coverage(bins, file)), .(cdition, file)]
  dat <- dcast(dat, 
               seqnames+start+end~cdition, 
               value.var = "score")
  saveRDS(dat, 
          "Rdata/gw_250bp_bins_quantif_cutnrun_replicates.rds")
}else if(!exists("dat"))
  dat <- readRDS("Rdata/gw_250bp_bins_quantif_cutnrun_replicates.rds")
dat <- na.omit(dat)

#-----------------------------------------------------------#
# Correlations
#-----------------------------------------------------------#
cols <- grep("rep", names(dat), value = T)
pdf("pdf/cutnrun/PCC_cutnrun_replicates_wo_ChIP-Seq.pdf")
vl_heatmap(cor(na.omit(dat[, ..cols])), 
           display_numbers = T, 
           legend_title = "PCC")
dev.off()

pdf("pdf/cutnrun/PCC_cutnrun_replicates_w_ChIP-Seq.pdf")
vl_heatmap(cor(na.omit(dat[, `H3K27Ac_ED_ChIP-Seq`:H3K27me3_PHD9_rep2])), 
           display_numbers = T, 
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

# with ChIP-Seq
pca <- as.data.table(prcomp(dat[, `H3K27Ac_ED_ChIP-Seq`:H3K27me3_PHD9_rep2])$rotation,
                     keep.rownames= T)
pca[grep(18, rn), Cc:= "grey50"]
pca[grep("D9", rn), Cc:= "gold"]
pca[grep("D11", rn), Cc:= "limegreen"]
pca[grep("29", rn), Cc:= "orangered"]
pca[grep("ChIP-Seq", rn), Cc:= "white"]
pca[grep("H3K27me3", rn), Cc:= colorRampPalette(c(Cc, "blue"))(3)[2], Cc]
pca[grep("H3K27Ac", rn), Cc:= colorRampPalette(c(Cc, "red"))(3)[2], Cc]

pdf("pdf/cutnrun/PCA_cutnrun_replicates_w_ChIP-Seq.pdf", 
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
                             normalizePath(list.files("db/bw/cutnrun_reps/", "^H3K27", full.names = T))))
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
                             normalizePath(list.files("db/bw/cutnrun_merge/", full.names = T))))
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

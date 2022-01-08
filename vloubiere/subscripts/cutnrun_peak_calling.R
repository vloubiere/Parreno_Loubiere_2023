setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
options(scipen = 999)
require(kohonen)
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)
require(BSgenome.Dmelanogaster.UCSC.dm6)

reads <- fread("Rdata/metadata_cutnrun_final.txt")
reads[, shuffled_bed:= paste0("db/bed/cutnrun/merge/", ChIP, "_", cdition, "_merge_shuffled.bed"), .(ChIP, cdition)]
reads[, peaks_narrowpeaks:= paste0("db/narrowpeaks/", ChIP, "_", cdition, "_peaks.narrowPeak"), .(ChIP, cdition)]
reads[, peaks_merged:= paste0("db/narrowpeaks/", ChIP, "_", cdition, "_merged_peaks.narrowPeak"), .(ChIP, cdition)]

# Peaks calling
reads[, {
  bs <- switch(ChIP, "H3K27Ac"= 100, "H3K27me3"= 1000)
  min_width <- switch(ChIP, "H3K27Ac"= 250, "H3K27me3"= 5000)
  min_dist <- switch(ChIP, "H3K27Ac"= 201, "H3K27me3"= 2501)
  min_qval <- switch(ChIP, "H3K27Ac"= 3, "H3K27me3"= 10)
  min_OR <- switch(ChIP, "H3K27Ac"= 1, "H3K27me3"= 2)
  # Create shuffled bed if not already done
  if(!file.exists(shuffled_bed))
    vl_exportBed(vl_shuffleBed(ChIP_bed_merge), shuffled_bed)
  # Peak calling
  if(!file.exists(peaks_narrowpeaks))
  {
    peaks <- vl_peakCalling(ChIP_bed_merge,
                            shuffled_bed,
                            gaussian_blur = T,
                            BSgenome = BSgenome.Dmelanogaster.UCSC.dm6, 
                            bins_width = bs, 
                            bins_OR_cutoff = 1)
    vl_exportBed(peaks, 
                 peaks_narrowpeaks)
  }
  # Merged peaks
  if(!file.exists(peaks_merged))
  {
    if(!exists(peaks))
      peaks <- vl_importBed(peaks_narrowpeaks)
    .m <- vl_collapseBed(peaks, mingap = min_dist)[end-start>min_width]
    merged <- vl_enrichBed(.m, 
                           ChIP_bed_merge, 
                           shuffled_bed)
    vl_exportBed(merged[qValue>min_qval & signalValue>=min_OR], 
                 peaks_merged)
  }
  print("done")
}, .(ChIP, cdition, ChIP_bed_merge, shuffled_bed, peaks_narrowpeaks, peaks_merged)]

# Segment genome into homogenous regions
if(!file.exists("Rdata/K27_regions_segmentation.bed"))
{
  merged_peaks <- reads[, vl_importBed(peaks_merged), .(ChIP, cdition, peaks_merged)]
  merged_peaks <- vl_collapseBed(merged_peaks, return_idx_only = T)
  setkeyv(merged_peaks, c("seqnames", "start"))
  K27 <- merged_peaks[, {
    ranges <- sort(unique(c(start, end[-.N]+1, end[.N])))
    .(start= ranges[-length(ranges)], end= ranges[-1]-1)
  }, .(seqnames, idx)]
  K27 <- unique(K27)[end-start>500]
  K27 <- cbind(K27,
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_shuffled.bed")[, .(K27Ac_WT= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PHD11_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed")[, .(K27Ac_PHD11= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PHD9_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed")[, .(K27Ac_PHD9= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27Ac_PH29_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed")[, .(K27Ac_PH29= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_shuffled.bed")[, .(K27me3_WT= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PHD11_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed")[, .(K27me3_PHD11= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PHD9_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed")[, .(K27me3_PHD9= signalValue)],
               vl_enrichBed(K27, 
                            ChIP_bed = "db/bed/cutnrun/merge/H3K27me3_PH29_merge_uniq.bed",
                            Input_bed = "db/bed/cutnrun/merge/H3K27me3_PH18_merge_uniq.bed")[, .(K27me3_PH29= signalValue)])
  fwrite(K27[, !"name"],
         "Rdata/K27_regions_segmentation.bed",
         sep= "\t", 
         col.names = T)
}

#-------------------#
# clustering
#-------------------#
# Import
dat <- fread("Rdata/K27_regions_segmentation.bed")
cols <- grep("^K27", names(dat))
dat[, (cols):= lapply(.SD, function(x) {
  lim <- quantile(x, c(0.001, 0.999))
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  log2(x)
}), .SDcols= cols]
blacklisted <- vl_importBed("/mnt/d/_R_data/genomes/dm6-blacklist.v2.bed")
dat[!(blacklisted[dat, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N>0)]
# Clip extremes
grid <- somgrid(xdim = 2, 
                ydim = 2, 
                topo = "hexagonal", 
                toroidal = T)
som <- supersom(list(as.matrix(dat[, .(K27Ac_WT, K27me3_WT)]),
                     as.matrix(dat[, .(K27Ac_PHD11, K27Ac_PHD9, K27Ac_PH29, K27me3_PHD11, K27me3_PHD9, K27me3_PH29)])), 
                grid = grid, 
                user.weights = c(10, 4))
              
dat[, cl:= som$unit.classif]
setkeyv(dat, "cl")

par(mfrow= c(1,2))
setcolorder(dat, 
            c("cl", "K27Ac_WT", "K27me3_WT"))
vl_heatmap(dat[, cl:K27me3_WT], 
           cluster_rows= F,
           cluster_cols= F, 
           breaks= c(-3,-0.25,0.25,3), 
           col= c("cornflowerblue", "white", "white", "tomato"))
setcolorder(dat, 
            c("cl", "K27Ac_PHD11", "K27Ac_PHD9", "K27Ac_PH29", "K27me3_PHD11", "K27me3_PHD9", "K27me3_PH29"))
vl_heatmap(dat[, cl:K27me3_PH29], 
           cluster_rows= F,
           cluster_cols= F, 
           auto_margins= F, 
           breaks= c(-3,-0.25,0.25,3), 
           col= c("cornflowerblue", "white", "white", "tomato"))




# coll <- vl_collapse_DT_ranges(enriched_bins,
#                               mingap = 300)
# coll
# enriched_bins
# peaks <- vl_bedEnrichment(coll,"db/bed/cutnrun/merge/H3K27Ac_PH18_merge_uniq.bed", BSgenome = BSgenome.Dmelanogaster.UCSC.dm6)

# pdf("test/screenshots_peak_calling.pdf", width = 50, height = 4)
# vl_screenshot(GRanges("chr2L:4928931-5062187"),
#               "db/bw/cutnrun_merge_vl/H3K27Ac_PH18_merge.bw", 
#               genome = "dm6",
#               highlight_regions = enriched_bins[padj<0.00001])
# dev.off()





# 
# 
# 
# 
# 
# # Filter_peaks
# filter <- peaks[, {
#   if(ChIP=="H3K27Ac")
#     width_cutoff <- 200
#   if(ChIP=="H3K27me3")
#     width_cutoff <- 1000
#   .SD[end-start>=width_cutoff & OR>2]
# }, .(ChIP, cdition)]
# 
# 
# if(!file.exists("Rdata/counts_cutnrun_100bp_bins_gw.txt"))
# {
#   # Count reads
#   bins <- vl_binBSgenome(BSgenome.Dmelanogaster.UCSC.dm6, bin_size = 100)
#   reads <- fread("Rdata/metadata_cutnrun_final.txt")
#   counts <- reads[, {
#     .c <- rbindlist(lapply(ChIP_bed, fread, select = 1:3, col.names = c("seqnames", "start", "end")))
#     res <- .c[bins, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
#     print("DONE")
#     cbind(bins, counts= res)
#   }, .(cdition= paste0(ChIP, "_", cdition))]
#   counts <- dcast(counts, 
#                   seqnames+start+end~cdition, 
#                   value.var = "counts")
#   fwrite(counts, "Rdata/counts_cutnrun_100bp_bins_gw.txt")
# }else if(!exists("counts"))
#   counts <- fread("Rdata/counts_cutnrun_100bp_bins_gw.txt")
# 
# # Identify regions enriched in at least one of the conditions
# regions <- melt(counts, c("seqnames", "start", "end"))
# regions[value>5, check:= scale(log2(value))>1.5, variable]
# regions[, check:= ifelse(any(check), T, F), .(seqnames, start, end)]
# regions <- dcast(regions[(check)],
#                  seqnames+start+end~variable, 
#                  value.var = "value")
# 
# # Compute FC pval
# dat <- melt(regions, 
#             id.vars = c("seqnames", "start", "end", "H3K27Ac_PH18", "H3K27me3_PH18"))
# dat[grepl("H3K27Ac", variable), 
#     c("ctl", "sum_value", "sum_ctl"):= .(H3K27Ac_PH18, sum(value), sum(H3K27Ac_PH18)), variable]
# dat[grepl("H3K27me3", variable), 
#     c("ctl", "sum_value", "sum_ctl"):= .(H3K27me3_PH18, sum(value), sum(H3K27me3_PH18)), variable]
# dat[, c("estimate", "pval"):= {
#   .f <- fisher.test(matrix(c(value, ctl, sum_value, sum_ctl), ncol= 2))
#   .f[c("estimate", "p.value")]
# }, .(value, ctl, sum_value, sum_ctl)]
# dat[, qval:= p.adjust(pval, "fdr")]
# dat[, log2OR:= log2(estimate)]
# dat[, filter:= any(qval<1e-5) & all(is.finite(log2OR)), seqnames:end]
# 
# # Clustering
# mat <- dcast(dat[(filter)], 
#              seqnames+start+end~variable, 
#              value.var = "log2OR")
# if(!exists("pdf/cutnrun/clustering_bins_wss.pdf"))
# {
#   set.seed(123)
#   k.max <- 15
#   wss <- sapply(1:k.max, 
#                 function(k){kmeans(mat[, H3K27Ac_PH29:H3K27me3_PHD9], 
#                                    k, 
#                                    nstart=50,
#                                    iter.max = 15 )$tot.withinss})
#   pdf("pdf/cutnrun/clustering_bins_wss.pdf")
#   plot(wss)
#   dev.off()
# }
# mat[, cl:= kmeans(mat[, H3K27Ac_PH29:H3K27me3_PHD9], 5)$cluster]
# 
# pdf("pdf/clustering/clustering_aggregate.pdf")
# agg <- mat[, c(lapply(.SD, mean), .N), cl, .SDcols= patterns("^H3")]
# agg[, cl:= paste0("cluster ", cl, " (", N, ")")]
# agg$N <- NULL
# vl_heatmap(agg)
# dev.off()
# 
# collapsed_regions <- mat[, vl_collapse_DT_ranges(copy(.SD), mingap = 301), cl, .SDcols= seqnames:end]
# 
# # Compute final FCs
# reads <- fread("Rdata/metadata_cutnrun_final.txt")
# reads[, total_reads:= sum(sapply(ChIP_bed, function(x) fread(cmd= paste0("wc -l ", x))$V1)), .(ChIP, cdition)]
# final <- reads[, {
#   .c <- rbindlist(lapply(ChIP_bed, fread, select = 1:3, col.names = c("seqnames", "start", "end")))
#   res <- .c[collapsed_regions, .N, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$N
#   print("DONE")
#   cbind(collapsed_regions, counts= res)
# }, .(cdition= paste0(ChIP, "_", cdition), total_reads)]
# final[, score:= counts/total_reads]
# FC <- dcast(final, 
#             seqnames+start+end+cl~cdition, 
#             value.var = "score")
# FC <- melt(FC, 
#            c("seqnames", "start", "end", "cl", "H3K27Ac_PH18", "H3K27me3_PH18"))
# FC[grep("K27Ac", variable), log2FC:= log2(value)-log2(H3K27Ac_PH18)]
# FC[grep("K27me3", variable), log2FC:= log2(value)-log2(H3K27me3_PH18)]
# FC <- dcast(FC, 
#             seqnames+start+end~variable, 
#             value.var = "log2FC")
# 
# if(!exists("pdf/cutnrun/clustering_final_wss.pdf"))
# {
#   set.seed(123)
#   k.max <- 15
#   wss <- sapply(1:k.max, 
#                 function(k){kmeans(FC[, H3K27Ac_PH29:H3K27me3_PHD9], 
#                                    k, 
#                                    nstart=50,
#                                    iter.max = 15 )$tot.withinss})
#   pdf("pdf/cutnrun/clustering_final_wss.pdf")
#   plot(wss)
#   dev.off()
# }
# FC[, cl:= kmeans(FC[, H3K27Ac_PH29:H3K27me3_PHD9], 5)$cluster]
# final_agg <- FC[, c(lapply(.SD, mean), .N), cl, .SDcols= patterns("^H3")]
# final_agg[, cl:= paste0("cluster ", cl, " (", N, ")")]
# final_agg$N <- NULL
# 
# vl_heatmap(final_agg)
# pheatmap::pheatmap(mat[, H3K27Ac_PH29:H3K27me3_PHD9], 
#                    kmeans_k = 6)
# mat[, coor:= paste0(seqnames, ":", start, "-", end)]
# mat <- as.matrix(mat[, H3K27Ac_PH29:H3K27me3_PHD9],
#                  mat$coor)
# mat <- mat[apply(apply(mat, 2, function(x) is.finite(x)), 1, all),]
# grid <-  somgrid(4, 4, "hexagonal", toroidal = T)
# set.seed(1)
# init <- mat[sample(nrow(mat), grid$xdim*grid$ydim),]
# som <- som(mat, 
#            grid, 
#            init= init)
# 
# # som.hc <- cutree(hclust(object.distances(som, "codes")), 8)
# som.hc <- kmeans(som$codes[[1]], centers = 5)$cluster
# 
# pdf("pdf/cutnrun/som_clustering.pdf", width = 9.5)
# par(mfrow= c(2, 3))
# for(name in colnames(mat))
# {
#   x <- som$codes[[1]][,name]
#   x[x>2] <- 2
#   x[x<(-2)] <- -2
#   Cc <- c("cornflowerblue", "white", "white", "tomato")
#   plot(som, 
#        type= "property", 
#        property= x, 
#        palette.name= colorRampPalette(Cc), 
#        zlim= c(-2,2), 
#        main= name, 
#        shape= "straight")
#   add.cluster.boundaries(som, som.hc)
# }
# dev.off()

# 
# ###############
# summary(som.wines)
# 
# ## xyf
# xyf.wines <- xyf(scale(wines), vintages, grid = somgrid(5, 5, "hexagonal"))
# summary(xyf.wines)
# 
# ## supersom example
# data(yeast)
# yeast.supersom <- supersom(yeast, somgrid(6, 6, "hexagonal"),
#                            whatmap = c("alpha", "cdc15", "cdc28", "elu"),
#                            maxNA.fraction = .5)
# 
# plot(yeast.supersom, "changes")
# 
# obj.classes <- as.integer(yeast$class)
# colors <- c("yellow", "green", "blue", "red", "orange")
# plot(yeast.supersom, type = "mapping", col = colors[obj.classes],
#      pch = obj.classes, main = "yeast data")
# [Package kohonen version 3.0.10 Index]
# 
# 
# 
# 
# stop()
# pdf("pdf/cutnrun/screenshots_peak_calling_QC.pdf")
# res[, {
#   set.seed(1)
#   vl_screenshot(GRanges("chr3R",
#                         IRanges(start = c(16.615e6, 6.625e6, 5e6), 
#                                 end = c(17.025e6, 7.1e6, 15e6))),
#                 list.files("db/bw/cutnrun_merge_vl/", full.names= T),
#                 highlight_regions = GRanges(.SD),
#                 genome = "dm6", 
#                 n_genes = 1)
#   print("")
# }, class]
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# cl <- dcast(dat, 
#             seqnames+start+end~variable, 
#             value.var = "log2OR")
# pheatmap::pheatmap(cl[, H3K27Ac_PH29:H3K27me3_PHD9], kmeans_k = 5)
# export(GRanges(res), "/mnt/c/Users/User/Desktop/test.bed")
# 
# pdf("pdf/cutnrun/screenshots_peak_calling_QC.pdf", width = 14)
# vl_screenshot(GRanges("chr3R", IRanges(13.8e6, 14.2e6)),
#               list.files("db/bw/cutnrun_merge/", "H3K27Ac", full.names= T),
#               highlight_regions = K27Ac,
#               genome = "dm6", 
#               n_genes = 1)
# vl_screenshot(GRanges("chr3R", IRanges(12.6e6, 17.04e6)),
#               list.files("db/bw/cutnrun_merge/", "H3K27me3", full.names= T),
#               highlight_regions = K27me3,
#               genome = "dm6", 
#               n_genes = 1)
# dev.off()
# 
# export(K27Ac, "db/peaks/K27Ac_cutnrun.bed")
# export(K27me3, "db/peaks/K27me3_cutnrun.bed")

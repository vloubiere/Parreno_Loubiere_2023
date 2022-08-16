setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)

# Import
FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
FBGN <- as.data.table(FBGN)
promoters <- FBGN[type=="gene" & gene_id %in% FBGN[type=="mRNA", gene_id], 
                  .(seqnames, 
                    start= ifelse(strand=="+", start, end), 
                    strand), .(gene_id, gene_symbol)]
promoters[, end:= start]

# Import peaks and K27 nonoverlapping peaks
peaks <- fread("Rdata/processed_metadata_CUTNRUN.txt")
peaks <- peaks[ChIP %in% c("H3K27me3", "H3K27Ac"), vl_importBed(peaks_merged), .(cdition, ChIP, peaks_merged)]
nonoverlapping <- vl_importBed("db/narrowpeaks/cutnrun/K27_segmentation_mergedBed.bed")
setkeyv(peaks, c("seqnames", "start", "end"))
setkeyv(nonoverlapping, c("seqnames", "start", "end"))

# Overlaps
peaks[, {
  .c <- foverlaps(.SD, nonoverlapping)
  .c <- .c[, .(ChIP= ChIP[which.max(min(c(end, i.end))-max(c(start, i.start)))]), .(seqnames, start, end)]
  nonoverlapping[.c, (cdition):= i.ChIP, on= c("seqnames", "start", "end")]
  print("")
}, cdition]
cols <- grep("^PH", names(nonoverlapping), value = T)
nonoverlapping[, (cols):= lapply(.SD, function(x) sapply(x, function(x) switch(x, "H3K27me3"= -1, "H3K27Ac"= 1, 0))), .SDcols= cols]
nonoverlapping <- na.omit(nonoverlapping[apply(nonoverlapping[, PHD9:PH29], 1, function(x) length(unique(x))>1)])

# Clustering
ChIP <- as.matrix(nonoverlapping[, K27me3_log2_enr_PH18:K27Ac_log2_enr_PH29])
ChIP_c <- apply(ChIP, 2, function(x){
  lim <- quantile(x, 0.95)
  x[x>lim] <- lim
  return(x)
})
layers <- list(enrich= ChIP_c,
               bin= as.matrix(nonoverlapping[, PHD9:PH29]))
grid <- somgrid(2,
                3,
                "hexagonal",
                toroidal = T)
N_cl <- grid$xdim*grid$ydim
set.seed(1)
som <- supersom(lapply(layers, scale), 
                grid= grid, 
                user.weights = c(1,2))
nonoverlapping[, cl:= som$unit.classif]

# Overlapping genes
cols <- grep("log2_enr", names(nonoverlapping), value = T)
nonoverlapping[, size:= apply(.SD, 1, function(x) diff(range(x))), .SDcols= cols]
genes <- promoters[nonoverlapping, .(cl= i.cl,
                                     size= i.size, 
                                     gene_id, 
                                     gene_symbol), .EACHI, on= c("seqnames", "start<=end", "end>=start")]
genes <- unique(na.omit(genes))

# Network
netGenes <- genes[cl %in% c(1,2,3,4,5)]
net <- vl_STRING_interaction(netGenes[, gene_symbol],
                             col = vl_palette_categ2(length(unique(netGenes$cl)))[netGenes[, cl]],
                             score_cutoff = 900,
                             size = netGenes[, log2(size+1)*2.5])

# Motif enrichment
vl_exportBed(nonoverlapping[cl==1, 1:3],
             "db/peaks/K27_cutnrun_changes/peaks_cluster_1_K27Ac_gain_transcient.bed")
vl_exportBed(nonoverlapping[cl==2, 1:3],
             "db/peaks/K27_cutnrun_changes/peaks_cluster_2_K27me3_loss_K27Ac_gain_PH29.bed")
vl_exportBed(nonoverlapping[cl==3, 1:3],
             "db/peaks/K27_cutnrun_changes/peaks_cluster_3_K27Ac_increase_PH29.bed")
vl_exportBed(nonoverlapping[cl==4, 1:3],
             "db/peaks/K27_cutnrun_changes/peaks_cluster_4_K27Ac_loss_PH29.bed")
vl_exportBed(nonoverlapping[cl==5, 1:3],
             "db/peaks/K27_cutnrun_changes/peaks_cluster_5_weak_K27me3_loss_weak_K27Ac_gain_all.bed")
vl_exportBed(nonoverlapping[cl==6, 1:3],
             "db/peaks/K27_cutnrun_changes/peaks_cluster_6_K27me3_loss_PH29.bed")

pdf("pdf/cutnrun/clustering_diff_regions_peak.pdf",
    width = 20, 
    height = 10)
layout(matrix(1:3, ncol= 3), widths = c(0.5,1,1))
par(mar= c(10,4,1,4))
vl_heatmap(layers[[1]][order(som$unit.classif),], 
           cluster_rows = F, 
           cluster_cols = F, 
           breaks = c(0, 1, 3, 5, 9), 
           col = c("black", "blue", "yellow", "orangered", "red"), 
           show_rownames = F, 
           auto_margins = F, 
           legend_title = "Enr. (log2)")

cl_vec <- length(som$unit.classif)-(cumsum(table(som$unit.classif))-table(som$unit.classif)/2)
axis(2, 
     at= cl_vec, 
     labels = names(cl_vec), 
     las= 1, 
     lwd= NA)
abline(h= length(som$unit.classif)-cumsum(table(som$unit.classif)), col= "white", lwd= 2)

# Overlapping genes
vl_GO_clusters(split(genes$gene_id, genes$cl), 
               N_top = 10, 
               cex.balloons = 2)
# Network
par(mar= c(2,2,2,2))
vl_STRING_network(net)
legend("topleft", 
       fill= vl_palette_categ2(length(unique(netGenes$cl))),
       legend = paste("cluster", seq(max(netGenes$cl))),
       bty= "n")
dev.off()

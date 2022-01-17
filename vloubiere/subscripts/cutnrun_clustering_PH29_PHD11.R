setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)

# Import gene 
if(!exists("FBGN"))
{
  FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
  GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
  FBGN <- as.data.table(FBGN)
  promoters <- FBGN[type=="gene", .(seqnames, start, end, strand), .(gene_id, gene_symbol)]
  promoters <- vl_resizeBed(promoters, center = "start", upstream = 5000, downstream = 2000, ignore.strand = F)
  setkeyv(FBGN, "gene_id")
}

# Import FC  data
K27 <- readRDS("Rdata/K27_regions_segmentation.rds")
K27[, ID:= paste0(seqnames, ":", start, "-", end)]
K27[, log2FoldChange:= log2(OR)]
K27[, cdition:= paste0(ChIP, "_", cdition)]


# Use only misregulated regions
dat <- K27[!grepl("PH18|PHD9", cdition), 
           .(check= any(padj<0.00001 & abs(log2FoldChange)>1), 
             cdition, 
             log2FoldChange, 
             padj), ID][(check)]
# dcast FC
mat <- dcast(dat, ID~cdition, value.var = "log2FoldChange")
mat <- na.omit(mat)
mat <- as.matrix(mat, 1)
# Clip outliers
mat[mat>3.5] <- 3.5
mat[mat<(-3.5)] <- -3.5
# dcast padj
pmat <- dcast(dat, ID~cdition, value.var = "padj")
pmat <- as.matrix(pmat, 1)
pmat <- pmat[match(rownames(mat), rownames(pmat)),]
# Use pmat to make bin table
bin <- pmat
bin[] <- 0
bin[pmat<0.05 & apply(mat, 2, between, 1, Inf)] <- 1 #sig up
bin[pmat<0.05 & apply(mat, 2, between, -Inf, -1)] <- -1 #sig down


# WT LVLS
WT_dat <- K27[ID %in% dat$ID & grepl("PH18", cdition), 
              .(cdition, 
                log2FoldChange, 
                padj), ID]
# dcast WT levels
WT <- dcast(WT_dat, ID~cdition, value.var = "log2FoldChange")
WT <- as.matrix(WT, 1)
WT <- WT[match(rownames(mat), rownames(WT)),]

# Clustering
seed <- 2
grid <- somgrid(10,10,"hexagonal",toroidal= T)
set.seed(seed)
som <- supersom(data= list(scale(WT),
                           bin,
                           scale(mat)), 
                grid= grid, 
                user.weights= c(2,5,5))
Nhc <- 5
hc <- cutree(hclust(dist(som$codes[[2]])), Nhc)[som$unit.classif]
# Interaction network
overlap_genes <- data.table(ID= rownames(mat), cl= hc)
overlap_genes <- cbind(overlap_genes, as.data.table(GenomicRanges::GRanges(rownames(mat))))
overlap_genes <- promoters[overlap_genes, .(ID, cl, gene_id, gene_symbol), 
                           .EACHI, on= c("seqnames", "end>=start", "start<=end")]
overlap_genes <- na.omit(overlap_genes)
overlap_genes[, size:= max(abs(mat[rownames(mat)==.BY[[1]],])), ID]
setorderv(overlap_genes, "size", order = -1)
overlap_genes <- overlap_genes[, .SD[1], gene_symbol]
Cc <- igraph::categorical_pal(Nhc)
net <- vl_STRING_interaction(symbols = overlap_genes$gene_symbol,
                             size = overlap_genes$size,
                             col= Cc[overlap_genes$cl],
                             score_cutoff = 900)
# SPlit diff expressed genes
diff <- dat[padj<0.05 & abs(log2FoldChange)>1]
diff <- split(diff$ID, paste0(diff$cdition, ifelse(diff$log2FoldChange>0, " UP", " DOWN")))
# SAVE
obj <- cbind(WT, mat, pmat)
colnames(obj) <- paste0(colnames(obj), c(rep("_log2FoldChange", 4), rep("_padj", 4)))
obj <- as.data.table(obj, keep.rownames = "coor", key= "coor")
obj[overlap_genes, overlapping_gene_symbols:= .(.(i.gene_symbol)), on= "coor==ID"]
obj[overlap_genes, overlapping_gene_FBgns:= .(.(i.gene_id)), on= "coor==ID"]
obj[, cl:= hc]
saveRDS(obj,
        "Rdata/clustering_cutnrun_PH29_PHD11.rds")

#------------------------------------------------#
# PLOT
#------------------------------------------------#
pdf("pdf/cutnrun/clustering_K27_cutnrun_PH29_PHD11.pdf", width = 21)
layout(matrix(c(1,2,3,4,5), nrow= 1), 
       widths = c(4,2,2,7,5))
# upset plot
vl_upset_plot(diff, 
              intersection_cutoff = 5)
# Clustering heatmap
vl_heatmap(WT[order(hc),],
           cluster_rows = F,
           cluster_cols = F,
           breaks = c(-1, 0, 4),
           col= c("cornflowerblue", "white", "red"), 
           show_rownames = F, 
           legend_title = "PH18 enrich. (log2)")
abline(h= nrow(mat)-cumsum(table(hc))+0.5)
# Clustering heatmap
vl_heatmap(mat[order(hc), c(1,3,2,4)],
           cluster_rows = F,
           cluster_cols = F,
           breaks = c(-4, -1, -0.25, 0.25, 1, 4),
           col= c("blue", "cornflowerblue", "white", "white", "tomato", "red"), 
           show_rownames = F, 
           legend_title = "log2FC over PH18")
abline(h= nrow(mat)-cumsum(table(hc))+0.5)
# GO
par(cex= 0.7)
vl_GO_clusters(FBgn_list = split(overlap_genes$gene_id, overlap_genes$cl),
               N_top = 10,
               cex.balloons = 1.5,
               padj_cutoff = 0.05)
# Network
par(mar= c(1,1,1,1), cex= 1)
vl_STRING_network(net,
                  cex.vertices = 4,
                  cex.vertices.labels = 0.5)
legend("topleft",
       bty= "n",
       fill= Cc,
       legend = paste0("cluster ", seq(max(hc))))
# Screenshots
layout(matrix(1))
par(mar= c(4,8,2,2))
obj[, order:= apply(.SD, 1, max), .SDcols= c("H3K27Ac_PH29_log2FoldChange", 
                                             "H3K27Ac_PHD11_log2FoldChange",
                                             "H3K27me3_PH29_log2FoldChange", 
                                             "H3K27me3_PHD11_log2FoldChange")]
setorderv(obj, "order", -1)
sel <- obj[, {
  set.seed(1)
  as.data.table(GenomicRanges::GRanges(.SD[1:3, coor]))
}, keyby= cl]
tmp <- tempfile(fileext = ".bed")
vl_exportBed(sel, tmp)
ext_sel <- sel[, .(seqnames, start= start-10000, end= end+10000)]
vl_screenshot_simple(bed = ext_sel,
                     space= 10,
                     gtf = "/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf",
                     tracks= c(tmp,
                               "db/bw/cutnrun_merge_vl/H3K27Ac_PH18_merge.bw",
                               "db/bw/cutnrun_merge_vl/H3K27Ac_PH29_merge.bw",
                               "db/bw/cutnrun_merge_vl/H3K27Ac_PHD11_merge.bw",
                               "db/bw/cutnrun_merge_vl/H3K27me3_PH18_merge.bw",
                               "db/bw/cutnrun_merge_vl/H3K27me3_PH29_merge.bw",
                               "db/bw/cutnrun_merge_vl/H3K27me3_PHD11_merge.bw"),
                     names= c("region",
                              "K27Ac 18",
                              "K27Ac 29",
                              "K27Ac D11",
                              "K27me3 18",
                              "K27me3 29",
                              "K27me3 D11"),
                     max = c(70,15,70,100,30,100))
x <- seq(par("usr")[1], par("usr")[2], length.out= nrow(sel)/3+1)
rect(x[-length(x)],
     par("usr")[3],
     x[-1],
     par("usr")[4],
     border= NA,
     col= adjustcolor(Cc, 0.2))
text(x[-length(x)]+diff(x)/2,
     par("usr")[4],
     paste0("cluster ", seq(Nhc)), 
     pos= 3)
dev.off()





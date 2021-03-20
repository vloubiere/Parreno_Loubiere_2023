dat <- readRDS("Rdata/som_clustering_transcriptomes.rds")
gtf <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
prom <- resize(gtf[gtf$gene_id %in% dat$FBgn & gtf$type=="gene"], 1, fix = "start")
prom <- resize(prom, 2000, "center")
seqlevelsStyle(prom) <- "UCSC"

sel <- BigWigSelection(prom, "score")
bw <- list.files("../../public_data/dm6/bw/", full.names = T)
bw <- bw[grepl("_ED", bw)]
bw <- bw[!grepl("INPUT", bw)]
cdition <- gsub("_merge|.bw", "" , basename(bw))
bw <- rbindlist(lapply(bw, function(x) as.data.table(import.bw(x, selection= sel))), idcol = "bw_file")
bw[, bw_file:= cdition[bw_file]]

.q <- dat[as.data.table(prom), , on= "FBgn==gene_id"]
.q <- unique(.q[, .(FBgn, cl, symbol, ycoor, seqnames, start, end , strand)])
.q <- .q[, .(bw_file= cdition), .q]
.q[, uniq_ID:= .I, .q]
res <- bw[.q, .(uniq_ID, bw_counts= mean(rep(score, width))), .EACHI, on= c("bw_file", "seqnames", "start<end", "end>start"), nomatch= 0]
.q[res, bw_counts:= i.bw_counts, on= "uniq_ID"]
.q[is.na(bw_counts), bw_counts:= 0]
.q[, norm_counts:= scale(log2(bw_counts+1)), bw_file]

setorderv(.q, "ycoor", -1)
my_heatmap(.q, row.BY = "ycoor", col.BY = "bw_file", value.var = "norm_counts", cluster_rows = F, breaks = c(-2,0, 2))
abline(h= .q[, max(ycoor), cl]$V1)

pdf("pdf/HTM_clusters.pdf", width = 8)
mat <- as.matrix(dcast(.q, cl~bw_file, value.var = "bw_counts", fun.aggregate = mean), 1)
my_heatmap(scale(mat), cluster_rows = F, leg_title = "mean zscore")
dev.off()

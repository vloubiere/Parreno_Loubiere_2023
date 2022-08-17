setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(cl)]
motifs <- fread("Rdata/final_RE_motifs_table.txt")
motifs <- dat[motifs, on="FBgn", nomatch= NULL]

# Compute gene network
sizes <- apply(dat[, ..cols], 1, function(x) max(abs(x), na.rm=T))
sizes[sizes>quantile(sizes, 0.99)] <- quantile(sizes, 0.99)
net <-vl_STRING_interaction(symbols = dat$symbol,
                      species = "Dm",
                      col= dat$col,
                      size = sizes*1.75,
                      cex.label = sizes/8,
                      plot= F)
# Compute GO enrich
GO_all <- vl_GO_enrich(split(dat$FBgn, 
                             dat$cl),
                        species = "Dm", 
                        plot= F)
GO_PRC1 <- vl_GO_enrich(split(dat$FBgn, 
                              dat[, .(cl, ifelse(PRC1_bound, "PRC1+", "PRC1-"))]),
                        species = "Dm", 
                        plot= F)

# Compute motif enrichments
groups <- split(motifs, motifs$group)
cols <- grep("_mot$", names(motifs), value = T)
enr <- lapply(groups, function(x) {
  vl_motif_cl_enrich(counts_matrix = as.matrix(x[, ..cols]),
                     cl_IDs = x[, paste0(cl, ifelse(PRC1_bound, " PRC1+", " PRC1-"))],
                     auto_margins= F, 
                     plot= F)
})

#################################
# PLOT
#################################
pdf("pdf/Figures/Clustering_GFP-_system_RNA.pdf", 
    width= 30)
mat <- matrix(1:7, nrow= 1, byrow = T)
layout(mat, 
       widths = c(1.5,2,2,2.25,1.5,1.5,1.5))

# Heatmap
par(mar= c(7,10,2,10))
cols <- c("log2FoldChange_PH29", "log2FoldChange_PHD9", "log2FoldChange_PHD11")
pl <- vl_heatmap(dat[, ..cols], 
                 row_clusters = dat$cl,
                 cluster_rows= F,
                 cluster_cols= F, 
                 breaks = seq(-5, 5, length.out= 100), 
                 show_rownames = F, 
                 col = vl_palette_blueWhiteRed(100), 
                 legend_title = "Fold Change (log2)", 
                 auto_margins = F,
                 show_col_clusters = F,
                 plot= F)
pl$rows[unique(dat[, .(cl, col)]), col:= i.col, on= "cl"]
plot(pl)
pl$rows[, {
  text(par("usr")[1], 
       mean(y),
       paste0(.N, " genes"),
       xpd= T,
       pos= 2)
}, cl]

# Network
par(mar= c(0,0,1,0),
    las= 0)
set.seed(1)
plot(net,
     score_cutoff= 600, 
     top_N= 1000)

leg <- unique(dat[, .(cl, col)])
legend("topleft",
       fill= leg$col,
       legend= paste0("Cluster ", leg$cl),
       bty= "n")
title("STRING interactions")

# GOs
par(las= 2,
    mar= c(4,20,1,12),
    mgp= c(3,0,0))
plot(GO_all,
     padj_cutoff = 0.05, 
     top_enrich = 5, 
     auto_margins= F,
     cex.balloons= 0.6)
title("GOs enrichment per cluster")
par(las= 2,
    mar= c(4,20,1,8),
    mgp= c(3,0,0))
plot(GO,
     padj_cutoff = 0.05, 
     top_enrich = 5, 
     auto_margins= F,
     cex.balloons= 0.4)
title("GOs enrichment per cluster +/- PRC1")

# Motifs
par(mar= c(4,6,1,8))
for(i in seq(enr))
{
  plot(enr[[i]],
       padj_cutoff= 0.05, 
       auto_margins= F, 
       top_enrich= c(5, 4, 5)[i], 
       cex.balloons= c(1.5, 0.5, 1)[i])
  title(main= names(groups)[i])
}

dev.off()


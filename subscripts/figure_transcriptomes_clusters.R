setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)


loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import data
obj <- readRDS("Rdata/clustering_RNA.rds")
list2env(obj, 
         envir = environment())

pdf("pdf/Figures/Clustering_GFP-_system_RNA.pdf", 
    width= 22)
mat <- matrix(c(1,2,5,6,7,
                1,3,5,6,7,
                1,4,5,6,7),
              nrow= 3, byrow = T)
layout(mat, 
       widths = c(1,0.9,1.5,2,1.3))

# Heatmap
par(mar= c(7,10,2,10))
cols <- grep("^log2FoldChange", names(data), value = T)
vl_heatmap(data[(order), ..cols], 
           row_clusters = data[(order), cl],
           cluster_rows= F,
           cluster_cols= F, 
           breaks = seq(-5, 5, length.out= 100), 
           show_rownames = F, 
           col = vl_palette_blueWhiteRed(100), 
           legend_title = "Fold Change (log2)", 
           auto_margins = F,
           show_col_clusters = F)
text(par("usr")[1], 
     cl_pos-diff(c(1, cl_pos))/2,
     rev(paste0("Cluster ", names(cl_counts), "\n(", cl_counts, ")")),
     xpd= T,
     pos= 2)
title("GFP- system (SOM)")

# CHROMHMM
all_chrom <- fread("external_data/chromatin_types_SA2020_table_s1.txt")
chrom <- rbind(data[, .(chromhmm, cl)],
               all_chrom[, .(chromhmm= V4, cl="all")])
chrom <- dcast(na.omit(chrom), 
               chromhmm~cl, 
               fun.aggregate = length)
chrom <- as.matrix(chrom, 1)
chrom <- apply(chrom, 2, function(x) x/sum(x)*100)
chrom <- chrom[order(apply(chrom, 1, sum), decreasing = T),]
Cc <- sapply(rownames(chrom), function(x) switch(x, 
                                                 "aTSS"= "tomato",
                                                 "aTTS"= "red",
                                                 "Enhancer"= "limegreen",
                                                 "PcG"= "cornflowerblue",
                                                 "Null"= "grey"))
par(mar= c(2,5,2,6))
barplot(chrom,
        col= Cc,
        las= 1,
        ylab= "% chromhmm class")
legend(par("usr")[2],
       par("usr")[4],
       legend = rownames(chrom),
       fill= Cc,
       xpd= T,
       bty= "n")

# PRC1 binding
all_PRC1 <- loadRData("external_data/SA2020_cl.list")
all_PRC1 <- rbindlist(lapply(all_PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
PRC1_binding <- rbind(data[, .(PRC1_cluster, cl)],
                      all_PRC1[, .(PRC1_cluster, cl= "all")])
PRC1_binding[, total:= ifelse(cl=="all", 16514, .N), cl]
PRC1_binding <- PRC1_binding[!is.na(PRC1_cluster)]
PRC1_binding[, PRC1_cluster:= switch(PRC1_cluster,
                                     "aTSS A"= "aTSS", 
                                     "aTSS B"= "aTSS", 
                                     "Enhancer A"= "Enhancer", 
                                     "Enhancer B"= "Enhancer", 
                                     "Polycomb A"= "Polycomb", 
                                     "Polycomb B"= "Polycomb", 
                                     "Repressed"= "Repressed"), PRC1_cluster]
PRC1_binding <- PRC1_binding[, .(perc= .N/total*100), .(cl, PRC1_cluster)]
mat <- dcast(unique(PRC1_binding), 
             PRC1_cluster~cl, 
             value.var = "perc", 
             fill = 0)
mat <- as.matrix(mat, 1)
mat <- mat[c("aTSS", "Enhancer", "Polycomb", "Repressed"),]
Cc <- c("tomato", "limegreen", "cornflowerblue", "grey")
bar <- barplot(mat,
               ylab= "% PRC1 bound genes",
               las= 1,
               col= Cc)
legend(par("usr")[2],
       par("usr")[4],
       legend = rownames(mat),
       fill= Cc,
       xpd= T,
       bty= "n")

# ref FPKMs
vl_boxplot(W18_FPKM~cl,
           data,
           las= 0,
           ylab= "W18_FPKM")

# Network
par(mar= c(0,0,1,0),
    las= 0)
plot(network, 
     score_cutoff= 600, top_N= 1000)
leg <- unique(data[network$V, .(cl, col), on= "symbol==name"])[order(cl)]
legend("topleft",
       fill= leg$col,
       legend= paste0("Cluster ", leg$cl),
       bty= "n")
title("STRING interactions")

# GOs
par(las= 2,
    mar= c(4,30,1,8),
    mgp= c(3,0,0))
plot(GO_cl_PRC1, 
     padj_cutoff= 0.05, 
     auto_margins= F, 
     top_enrich = 5, 
     cex.balloons= 0.5)
title("GOs enrichment")

# Motifs
par(mar= c(4,10,1,8))
plot(mot_enrichment, 
     padj_cutoff= 0.05, 
     auto_margins= F, 
     top_enrich= 4,
     cex.balloons= 1.5)
title("Motifs enrichment")
dev.off()


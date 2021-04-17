# Clustering heatmaps
dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")

pdf("pdf/clustering_transcriptomes.pdf", width = 10, height = 10)
par(mfrow= c(1,2), mar= c(8.1,8.1,5.1,6.1))

Cc <- c("blue", "cornflowerblue", "white", "tomato", "red")
mat <- as.matrix(dcast(dat, 
                       cl+FBgn~cdition, 
                       value.var = "clipped_log2FoldChange")[,-2], 1)
par(xaxs= "i", 
    yaxs= "i", 
    las= 2)
vl_heatmap(mat, 
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "log2FC")
box()
cllines <- 1-cumsum(dat[, .N, keyby= cl][,N])/nrow(dat)
abline(h= cllines[-length(cllines)])
axis(2, 
     labels = dat[, cl, keyby= cl][, cl],
     at = cllines-diff(c(1, cllines))/2, 
     tick = 0)


# Plot aggregate
agg <- as.matrix(dcast(dat, 
                       cl~cdition, 
                       fun.aggregate = mean,
                       value.var = "clipped_log2FoldChange"), 1)
vl_heatmap(agg, 
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "log2FC")
box()
axis(2, 
     labels = dat[, paste0("cluster ", cl, " (", length(unique(FBgn)), ")"), keyby= cl][, V1],
     at = 1-dat[, cl, keyby= cl][, cl]/length(unique(dat$cl))+1/(2*length(unique(dat$cl))), 
     tick = 0)

dev.off()


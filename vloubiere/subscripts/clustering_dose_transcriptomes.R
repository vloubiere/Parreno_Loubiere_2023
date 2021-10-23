setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
dir.create("pdf/clustering", 
           showWarnings = F)


dat <- readRDS("Rdata/final_FC_table.rds")
mat <- dcast(dat[cdition %in% c("PH18_dose", "PH21_dose", "PH25_dose", "PH29_dose", "PH29_TS")], 
             symbol~cdition, 
             value.var = "log2FoldChange")
mat <- as.matrix(mat, 1)
mat <- na.omit(mat)
mat <- mat[rownames(mat) %in% dat[padj<0.05 & abs(log2FoldChange)>1, symbol],]
# Clip outliers
mat <- apply(mat, 2, function(x)
{
  lim <- quantile(x, c(0.01, 0.99), na.rm= T)
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  return(x)
})

# Elbow Method for finding the optimal number of clusters
set.seed(123)
k.max <- 15
wss <- sapply(1:k.max, 
              function(k){kmeans(mat, 
                                 k, 
                                 nstart=50,
                                 iter.max = 15 )$tot.withinss})

# PLOT
pdf("pdf/clustering/clustering_dose_transcriptomes.pdf", 
    width = 4,
    height = 4.5)
plot(wss, 
     xlab= " Number of clusters",
     ylab= "wss")
lines(wss)
points(5, 
       wss[5], 
       cex= 2, 
       col= "red")
par(mai= c(1.5651667, 0.75000000, 0.4500000, 0.9371667))
cl <- vl_heatmap(mat, 
                 cluster_cols = T, 
                 show_rownames = F, 
                 breaks = c(-5,0,5),
                 col =  c("darkblue", "blue", "cornflowerblue", "white", "tomato", "red", "darkred"),
                 kmeans_k = 5, 
                 legend_title = "log2FC", 
                 auto_margins = F)
cl[, text(0, 
          mean(range(y))/max(cl$y), 
          pos= 2,
          offset= 0.25,
          paste0(formatC(length(unique(row)), big.mark = ","), " genes"), 
          xpd= T, 
          cex= 0.6), rcl]
rect(0,0,1,1,xpd=T)
dev.off()

# SAVE
cl[, symbol:= row]
cl[unique(dat[, .(symbol, FBgn)]), FBgn:= i.FBgn, on= "symbol"]
saveRDS(cl, "Rdata/clustering_dose_transcriptomes.rds")

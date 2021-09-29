setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
dir.create("pdf/clustering", 
           showWarnings = F)
#-------------------------------#
# Import cl data and recreate matrix for plotting
#-------------------------------#
dat <- readRDS("Rdata/clustering_dose_transcriptomes.rds")
cl_mat <- dcast(dat, row~col, value.var = "value")
cl_mat <- as.matrix(cl_mat, 1)

#-------------------------------#
# Compute features per gene and make matrix
#-------------------------------#
coll <- unique(dat[, .(rcl, FBgn, y)])
setorderv(coll, "y", -1)
feat <- readRDS("Rdata/ED_REs_features_final.rds")

res <- feat[coll, .(y,
                    PH= max(c(PH_score), na.rm= T),
                    H2AK118Ub= max(c(H2AK118Ub_score), na.rm= T),
                    SUZ12= max(c(SUZ12_score), na.rm= T),
                    H3K27me3= max(c(H3K27me3_score), na.rm= T),
                    H3K27ac= max(c(H3K27ac_score), na.rm= T),
                    H3K27me2= max(c(H3K27me2_score), na.rm= T)), .EACHI, on= "FBgn"]
mat <- as.matrix(res[,PH:H3K27me2])
mat[mat==(-Inf)] <- NA
mat <- apply(mat, 2, scale)

pdf("pdf/clustering/features_dose_clusters.pdf", 
    width = 8, 
    height = 4.5)
par(mfrow=c(1, 2),
    mai= c(1.5651667, 0.75000000, 0.4500000, 1.2))
vl_heatmap(cl_mat, 
           cluster_cols = T, 
           show_rownames = F, 
           breaks = c(-5,0,5),
           col =  c("darkblue", "blue", "cornflowerblue", "white", "tomato", "red", "darkred"),
           kmeans_k = 5, 
           legend_title = "log2FC", 
           auto_margins = F)
vl_heatmap(mat, 
           cluster_rows = F, 
           cluster_cols = F, 
           breaks = c(-2, 0, 5), 
           auto_margins = F, 
           legend_title = "ChIP enrichment", 
           show_rownames = F)
.l <- dat[, max(y)/nrow(cl_mat), rcl][, V1]
segments(x0 = 0,
         y0 = .l,
         x1= 1,
         y1= .l)
dev.off()

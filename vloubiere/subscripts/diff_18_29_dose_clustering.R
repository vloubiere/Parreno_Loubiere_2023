setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
dir.create("pdf/comparison_ph18_ph29_DOSE", 
           showWarnings = F)
require(igraph)
require(STRINGdb)

#-------------------------#
# Import ph29 vs ph18 dose comparison
#-------------------------#
dat <- data.table(file= c("db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH18_vs_RNA_W18.txt",
                          "db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH29_vs_RNA_W29.txt"))
dat[, cdition:= c("PH18", "PH29")]
dat <- dat[, fread(file), (dat)]
names(dat)[names(dat)=="V1"] <- "FBgn"
# Add Symbols
symbols <- unique(readRDS("Rdata/final_FC_table.rds")[, 2:3])
dat[symbols, symbol:= i.symbol, on= "FBgn"]
# Select genes which are differentially affected in the two conditions
sel <- dat[, (length(which(padj<0.05))==1 | (all(padj<0.05) & sign(log2FoldChange[1])!=sign(log2FoldChange[2]))) 
           & any(abs(log2FoldChange)>1), FBgn][(V1), FBgn]
sub <- dat[FBgn %in% sel]

#-------------------------#
# make plotting objects
#-------------------------#
# Heatmap
mat <- dcast(sub, FBgn~cdition, value.var = "log2FoldChange")
mat <- as.matrix(mat, 1)

# clusters
cl <- vl_heatmap(mat, 
                 kmeans_k = 4, 
                 breaks = c(-5,0,5),
                 plot= F)
cl <- cl[, .(size= abs(mean(value))), .(FBgn= row, cl= rcl)]
cl[, Cc:= c("tomato", "gold", "limegreen", "cornflowerblue")[.GRP], cl]
cl[symbols, symbol:= i.symbol, on= "FBgn"]

# network clusters legend

# network
net <- vl_STRING_interaction(cl$symbol, 
                             size = cl$size, 
                             cex.label = cl$size, 
                             col = cl$Cc,
                             score_cutoff = 500)
# FBgn list GOs
.l <- split(cl$FBgn, cl$cl)


#-------------------------#
# PLOT
#-------------------------#

pdf("pdf/comparison_ph18_ph29_DOSE/all_plots.pdf", 
    width = 16, 
    height = 7)
layout(matrix(1:3, ncol= 3), 
       widths = c(3,1,3))

vl_STRING_network(net, 
                  cex.vertices = 2,
                  cex.vertices.labels = 0.4, 
                  vertex.border.col = NA) 

legend("topleft", 
       pch= 19,
       col= unique(cl$Cc),
       legend = paste0("cl ", unique(cl$cl)), 
       bty= "n",
       cex= 2)

par(mar=c(2.5,3,7,3))
res <- vl_heatmap(mat,
                  kmeans_k = 4, 
                  breaks = c(-5,0,5), 
                  show_rownames = F, 
                  auto_margins = F, 
                  legend_title = "log2FC")
res[, text(0, 
          mean(range(y))/max(res$y), 
          pos= 2,
          offset= 0.25,
          paste0(formatC(length(unique(row)), big.mark = ","), " genes"), 
          xpd= T, 
          cex= 0.6), rcl]

vl_GO_clusters(.l, 
               N_top = 10, 
               padj_cutoff = 0.05, 
               cex = 0.6)
dev.off()

saveRDS(res, 
        "Rdata/clustering_dose_18_29_transcriptomes.rds")


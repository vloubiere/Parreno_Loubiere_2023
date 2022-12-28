# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
dat <- fread("Rdata/final_gene_features_table.txt")

genes <- list("JAK-STAT"= c("Socs36E", "Stat92E", "chinmo", "upd1", "upd2", "upd3", "zfh1"),
              "JNK"= c("Atf3", "Irbp18", "Ets21C", "Jra", "Mmp1", "Pax", "Pdp1", "cher", "ftz-f1", "kay", "puc"))
genes <- rbindlist(lapply(genes, function(x) data.table(symbol= x)), 
                   idcol = "group")
genes <- dat[genes, on="symbol"]
mat <- genes[, .(symbol, log2FoldChange_PH18, log2FoldChange_PH29, log2FoldChange_PHD9, log2FoldChange_PHD11)]
mat <- as.matrix(mat, 1)

pdf("pdf/Figure_2_heatmap_FC_JAKSTAT_JNK.pdf", 
    width = 3.5,
    height = 4.5)
par(mar= c(10,7,1,6),
    las= 2, 
    cex.axis= 0.8)
vl_heatmap(mat, 
           row_clusters = factor(genes$group, unique(genes$group)),
           cluster_rows= F,
           cluster_cols= F, 
           display_numbers= T,
           legend_title= "log2FC",
           display_numbers_matrix = round(mat, 1),
           display_numbers_cex = 0.6,
           breaks= c(-4, 0, 4), 
           show_colnames = F)
vl_tilt_xaxis(1:4, 
              labels = c("no ph-KD", "Constant ph-KD", "Transient ph-KD d9", "Transient ph-KD d11"))
dev.off()
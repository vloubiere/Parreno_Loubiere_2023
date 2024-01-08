setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[, .(symbol, log2FoldChange_PH18, log2FoldChange_PH29, log2FoldChange_PHD9, log2FoldChange_PHD11)]

genes <- data.table(symbol= c("upd1", "upd2", "upd3", "dome", "hop", "Stat92E", "zfh1", "chinmo", "Socs36E"),
                    group= c("Ligand", "Ligand", "Ligand", "Receptor", "Receptor TyrK", "TF", "Target genes", "Target genes", "Target genes"))
genes[, group:= factor(group, unique(group))]
genes <- dat[genes, on= "symbol"]
mat <- as.matrix(genes[, symbol:log2FoldChange_PHD11], 1)

# plot ----
pdf("pdf/Figure_2_heatmap_FC_JAKSTAT_JNK.pdf", 
    width = 3,
    height = 5)
par(mai= c(2,2,1,.2),
    las= 2, 
    cex.axis= 11/12)
vl_heatmap(mat, 
           cluster.rows= F,
           cluster.cols= F, 
           display.numbers= T,
           legend.title= "log2FC",
           display.numbers.matrix = round(mat, 1),
           display.numbers.cex = 0.5,
           breaks= c(-4, 0, 4), 
           show.rownames = T,
           show.colnames = F,
           tilt.colnames = T)
axis(1, 1:4, c("No ph-KD", "Constant ph-KD", "Transient ph-KD d9", "Transient ph-KD d11"), lwd= 0)
dev.off()
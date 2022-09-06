setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_noGFP"]
meta <- na.omit(meta[, .(cdition, FC_file)])
meta[, cdition:= factor(cdition, 
                        levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dat <- meta[, fread(FC_file), (meta)]
genes <- fread("Rdata/list_genes_interest_heatmaps.txt")
pl <- dat[genes, on= "FBgn", nomatch= NULL]
pl <- pl[!is.na(cdition)]
pl <- dcast(pl, 
            group+symbol+FBgn~cdition, 
            value.var = "log2FoldChange")
setorderv(pl, c("group", "symbol", "FBgn"))
mat <- as.matrix(pl[, .(symbol, PH18, PH29, PHD9, PHD11)], 1)

pdf("pdf/Figures/Heatmap_FC_GOF.pdf", 
    height = 50,
    width = 4)
par(mar= c(7,7,3,6),
    las= 1)
.c <- vl_heatmap(mat, 
                 row_clusters = pl$group,
                 cluster_rows= F,
                 cluster_cols= F, 
                 display_numbers= T,
                 legend_title= "log2FC",
                 display_numbers_matrix = round(mat, 1),
                 display_numbers_cex = 0.7,
                 breaks= c(-5, 0, 5))
dev.off()

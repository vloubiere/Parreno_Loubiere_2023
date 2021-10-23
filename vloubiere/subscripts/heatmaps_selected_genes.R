setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

dat <- readRDS("Rdata/final_FC_table.rds")
genes <- as.data.table(read_excel("Rdata/list_genes_interest.xlsx"))
genes[, group:= factor(group, levels= c("PRC1",
                                        "PRC2",
                                        "HOX",
                                        "JAK-STAT",
                                        "Cancer markers",
                                        "RDGN",
                                        "Ecdysone metamorhphosis",
                                        "cell-cell adhesion GO"))]
setorderv(genes, "group")

# epiCancer
sub <- dat[grepl("TS", cdition)]
sub <- dcast(sub, 
             symbol~cdition, 
             value.var = "log2FoldChange")
setkeyv(sub, "symbol")

pdf("pdf/hypothesis_driven_heatmaps_epiCancer.pdf", 
    height = 16*0.8, 
    width = 9*0.8)
par(xaxs= "i",
    yaxs= "i")
layout(matrix(c(1,2,3,4,5,6,7,8,8,8,8,8,8,8), 
              ncol= 2), 
       heights = c(6/36, 5/36, 11/36, 4/36, 5/36, 6/36, 9/36), 
       widths = c(0.8,1))
genes[, {
  mat <- as.matrix(na.omit(sub[symbol]), 1)
  par(mai= c(ifelse(group %in% c("Ecdysone metamorhphosis",
                                 "cell-cell adhesion GO"), 0.6, 0.1),
             0.75,
             0.3,
             ifelse(group=="cell-cell adhesion GO", 1, 0.15)))
  vl_heatmap(mat,  
             cluster_cols = F,
             cluster_rows = T,
             show_colnames = ifelse(group %in% c("Ecdysone metamorhphosis",
                                                 "cell-cell adhesion GO"), T, F),
             legend_title = "log2FC",
             auto_margins = F, 
             breaks = c(-5,0,5), 
             display_numbers = T)
  box(lwd= 0.25)
  title(group, line = 1)
}, group]
dev.off()

# Dose
sub <- dat[grepl("dose", cdition)]
sub <- dcast(sub, 
             symbol~cdition, 
             value.var = "log2FoldChange")
setkeyv(sub, "symbol")

pdf("pdf/hypothesis_driven_heatmaps_dose.pdf", 
    height = 16*0.8, 
    width = 9*0.8)
par(xaxs= "i",
    yaxs= "i")
layout(matrix(c(1,2,3,4,5,6,7,8,8,8,8,8,8,8), 
              ncol= 2), 
       heights = c(6/36, 5/36, 11/36, 4/36, 5/36, 6/36, 9/36), 
       widths = c(0.8,1))
genes[, {
  mat <- as.matrix(na.omit(sub[symbol]), 1)
  par(mai= c(ifelse(group %in% c("Ecdysone metamorhphosis",
                                 "cell-cell adhesion GO"), 0.6, 0.1),
             0.75,
             0.3,
             ifelse(group=="cell-cell adhesion GO", 1, 0.15)))
  vl_heatmap(mat,  
             cluster_cols = F,
             cluster_rows = T,
             show_colnames = ifelse(group %in% c("Ecdysone metamorhphosis",
                                                 "cell-cell adhesion GO"), T, F),
             legend_title = "log2FC",
             auto_margins = F, 
             breaks = c(-5,0,5), 
             display_numbers = T)
  box(lwd= 0.25)
  title(group, line = 1)
}, group]
dev.off()

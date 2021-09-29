setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

dat <- readRDS("Rdata/final_FC_table.rds")
setkeyv(dat, "symbol")
genes <- read_excel("Rdata/list_genes_interest.xlsx")

pdf("pdf/hypothesis_driven_heatmaps.pdf", 
    height = 4, 
    width = 6)
par(mar= c(5,7,4,5))

for(i in seq(genes))
{
  .c <- genes[[i]]
  mat <- dcast(na.omit(dat[.c]), 
               symbol~cdition,
               value.var = "log2FoldChange")
  mat <- as.matrix(na.omit(mat), 1)
  mat <- mat[na.omit(match(.c, rownames(mat))),]
      
  pl <- vl_heatmap(mat,
                   display_numbers = T,
                   cluster_cols = F,
                   cluster_rows = T, 
                   breaks = c(-1, 0, 2))
  pl[dat, padj:= i.padj, on= c("row==symbol", "col==cdition")]
  pl[, padj:= cut(padj, 
                  c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf),
                  c("****", "***", "**", "*", ""))]
  text(x = pl$xplot,
       y = pl$yplot, 
       labels = pl$padj, 
       pos= 3, 
       offset = 0.25,
       cex= 0.5)
}
dev.off()

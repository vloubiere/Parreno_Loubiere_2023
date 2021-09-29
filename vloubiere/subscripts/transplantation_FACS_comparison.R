library(ComplexUpset)
require("PerformanceAnalytics")

dat <- data.table(file= list.files("db/FC_tables/", full.names = T))
dat <- dat[grepl("_T_", file)]
dat <- dat[, fread(file), file]
dat[, cdition:= gsub("RNA_epiCancer_|_FC.txt", "", basename(file)), file]

pdf("pdf/upset_plot_comparison_up_genes_transplant.pdf")
intersect <- split(dat[grepl("vs_PH29$|vs_PHJ11$", cdition) & padj<0.001 & log2FoldChange>0, .(cdition, V1)], by= "cdition")
intersect <- lapply(intersect, function(x) x$V1)
vl_upset_plot(intersect)
dev.off()

pdf("pdf/upset_plot_comparison_down_genes_transplant.pdf")
intersect <- split(dat[grepl("vs_PH29$|vs_PHJ11$", cdition) & padj<0.001 & log2FoldChange<0, .(cdition, V1)], by= "cdition")
intersect <- lapply(intersect, function(x) x$V1)
vl_upset_plot(intersect)
dev.off()


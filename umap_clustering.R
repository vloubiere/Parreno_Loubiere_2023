setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(data.table)
require(umap)
require(circlize)

dat <- data.table(file= list.files("/groups/stark/vloubiere/projects/epigenetic_cancer/FC_tables", full.names = T))
dat <- dat[grep("PH18_vs_W18|PHJ9_vs_WKD|PHJ11_vs_WKD|PH29_vs_W29|PHJ9_vs_PH29|PHJ11_vs_PH29", file)]
dat[, cdition := gsub("_FC.txt", "", basename(file))]
dat <- dat[, fread(file), c(colnames(dat))]

dmat <- dcast(dat, rn~cdition, value.var = "log2FoldChange")
dmat[dat[cdition==cdition[1]], baseMean:= log2(i.baseMean+1), on= "rn"]
check <- dat[, any(padj<0.01), rn][(V1), rn]
dmat <- dmat[rn %in% check]

mat <- as.matrix(dmat, 1)
FCumap <- umap(mat)
Cc <- colorRamp2(c(-3, 0, 3), c("blue", "lightgrey", "red"))
plot(FCumap$layout, col= adjustcolor(Cc(mat$PH29_vs_W29), 0.5), pch= 16, cex= 0.5)

plot(FCumap$layout, col= adjustcolor(Cc(dmat$baseMean), 0.5), pch= 16, cex= 0.5)


plot(FCumap$layout, col= adjustcolor(Cc(dmat$PH18_vs_W18), 0.5), pch= 16, cex= 0.5)
plot(FCumap$layout, col= adjustcolor(Cc(dmat$PHJ11_vs_PH29), 0.5), pch= 16, cex= 0.5)
plot(FCumap$layout, col= adjustcolor(Cc(dmat$PHJ9_vs_WKD), 0.5), pch= 16, cex= 0.5)
plot(FCumap$layout, col= matlab.like2(15)[kmeans(FCumap$knn$distances, 15)$cluster], pch= 16, cex= 0.5)

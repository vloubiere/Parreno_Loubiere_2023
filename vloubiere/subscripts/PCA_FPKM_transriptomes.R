setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
dat <- lapply(c("db/dds/RNA/epiCancer_ED_allograft_RNA_gDNA_dds.rds",
                "db/dds/RNA/epiCancer_ED_RNA_CUTNRUN_dds.rds"), function(x) 
         {
           dat <- DESeq2::fpkm(readRDS(x))
           dat <- melt(as.data.table(dat, keep.rownames = T), id.vars = "rn")
           dat[, variable:= gsub("_1.bam|_2.bam|_3.bam", "", variable)]
           dat <- dat[, .(FPKM= mean(value)), .(rn, variable)]
         })
names(dat) <- c("GFP+_system", "GFP-_system")
dat <- rbindlist(dat, idcol = T)
dat[, variable:= gsub("^TPH", "PH", variable)]
dat[, variable:= gsub("^TW", "W", variable)]
dat[variable=="PHJ9", variable:= "PHD9"]
dat[variable=="PHJ11", variable:= "PHD11"]
dat[, col:= vl_palette_many_categ(12)[.GRP], variable]
dat <- dat[, .id:= paste0(variable, "_", .id)]
Cc <- unique(dat[, .(.id, col)])
dat <- dcast(dat, rn~.id, value.var = "FPKM")
mat <- as.matrix(dat, 1)
mat[is.na(mat)] <- 0
mat <- scale(mat)
pca <- prcomp(mat)$rotation

pdf("pdf/Figures/PCA_FPKMs_RNA.pdf",
    width = 5,
    height = 5.5)
plot(pca[,"PC1"],
     pca[,"PC2"],
     col= adjustcolor(Cc[rownames(pca), col, on= ".id"], 0.7),
     pch= ifelse(grepl("GFP\\+", rownames(pca)), 15, 19),
     las= 1,
     cex= 1.5)
legend("topright", 
       rownames(pca),
       pch= ifelse(grepl("GFP\\+", rownames(pca)), 15, 19),
       col= adjustcolor(Cc[rownames(pca), col, on= ".id"], 0.7),
       bty= "n",
       cex= 0.7)
dev.off()



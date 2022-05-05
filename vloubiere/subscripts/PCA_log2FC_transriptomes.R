setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[project=="RNA_epiCancer"]
meta <- meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), .(DESeq2_object, cdition)]
meta <- unique(meta[FC_file!="NA" & !grepl("vs_RNA_PHD11.txt$", FC_file)]) # Do not consider PHD11_T vs PHD11
meta <- meta[, fread(FC_file), (meta)]
meta[, cdition:= paste(ifelse(grepl("allograft", DESeq2_object), "allograft", "cutnrun"), cdition)]
meta[, cdition:= factor(cdition, levels= c("allograft RNA_PH18",
                                           "allograft RNA_PHD9",
                                           "allograft RNA_PHD11",
                                           "allograft RNA_PH29",
                                           "allograft RNA_PHD11_T2",
                                           "allograft RNA_PHD11_T3",
                                           "cutnrun RNA_PH18",
                                           "cutnrun RNA_PHD9",
                                           "cutnrun RNA_PHD11",
                                           "cutnrun RNA_PH29"))]
dat <- dcast(meta, FBgn~cdition, value.var = "log2FoldChange")
dat <- na.omit(dat)
pca <- as.data.table(prcomp(as.matrix(dat, 1))$rotation, 
                     keep.rownames = "cdition")
pca[grepl("allograft", cdition), col:= vl_palette_many_categ(10)[.GRP], cdition]
pca[grepl("cutnrun", cdition), col:= vl_palette_many_categ(10)[.GRP], cdition]


pdf("pdf/Figures/PCA_log2FC_RNA.pdf", 
    width = 5, 
    height = 5.5)
par(las= 1)
plot(pca$PC1,
     pca$PC2, 
     bg= adjustcolor(pca$col, 0.5),
     pch= ifelse(grepl("allograft", pca$cdition), 21, 22),
     cex= 2,
     las= 1, 
     xlab= "PC1",
     ylab= "PC2")
legend("topright",
       pt.bg= adjustcolor(pca$col, 0.5),
       legend= pca$cdition,
       bty= "n",
       pch= ifelse(grepl("allograft", pca$cdition), 21, 22))
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[project=="RNA_epiCancer"]
meta <- meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), .(DESeq2_object, cdition)]
meta <- unique(meta[FC_file!="NA"])
meta <- meta[, fread(FC_file), (meta)]
meta[, cdition:= paste(fcase(grepl("GFP\\+", DESeq2_object), "GFP+",
                             grepl("GFP\\-", DESeq2_object), "GFP-"), cdition)]
meta[, cdition:= factor(cdition, levels= c("GFP+ RNA_PH18",
                                           "GFP+ RNA_PHD9",
                                           "GFP+ RNA_PHD11",
                                           "GFP+ RNA_PH29",
                                           "GFP- RNA_PH18",
                                           "GFP- RNA_PHD9",
                                           "GFP- RNA_PHD11",
                                           "GFP- RNA_PH29"))]
dat <- dcast(meta, 
             FBgn~cdition, 
             value.var = "log2FoldChange")
dat <- na.omit(dat)
pca <- as.data.table(prcomp(as.matrix(dat, 1))$rotation, 
                     keep.rownames = "cdition")
pca[grepl("GFP+", cdition), col:= vl_palette_many_categ(4)[.GRP], cdition]
pca[grepl("GFP-", cdition), col:= vl_palette_many_categ(4)[.GRP], cdition]


pdf("pdf/Figures/PCA_log2FC_RNA.pdf", 
    width = 5, 
    height = 5.5)
par(las= 1)
plot(pca$PC1,
     pca$PC2, 
     bg= adjustcolor(pca$col, 0.5),
     pch= ifelse(grepl("GFP-", pca$cdition), 21, 22),
     cex= 2,
     las= 1, 
     xlab= "PC1",
     ylab= "PC2")
legend("bottomright",
       pt.bg= adjustcolor(pca$col, 0.5),
       legend= pca$cdition,
       bty= "n",
       pch= ifelse(grepl("allograft", pca$cdition), 21, 22))
dev.off()

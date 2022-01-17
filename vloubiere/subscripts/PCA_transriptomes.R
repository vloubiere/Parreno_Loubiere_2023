setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(diagram)

dat <- rbindlist(readRDS("Rdata/RNA_tables_object.rds")$FC, idcol = T)
dat[, cdition:= paste0(cdition, "_", .id)]
dat <- dcast(dat,
             FBgn~cdition, 
             value.var = "log2FoldChange")
dat <- na.omit(dat)


pca <- as.data.table(prcomp(as.matrix(dat, 1))$rotation, 
                     keep.rownames = "cdition")
Cc <- c("orchid1", "darkorchid1", "purple", "darkorchid4", 
        "olivedrab1", "limegreen", "olivedrab3", "olivedrab4", 
        "lightsteelblue1", "cornflowerblue", "blue", "navy", 
        "pink", "pink3", "indianred2", "red2", 
        "black", "gold", "goldenrod", "goldenrod4", "sienna4")

pdf("pdf/RNA_timecourse/PCA_transcriptomes.pdf", width = 5, height = 5.5)
par(las= 1)
plot(pca$PC1,
     pca$PC2, 
     bg= adjustcolor(colorRampPalette(Cc)(nrow(pca)), 0.5),
     pch= ifelse(grepl("allograft", pca$cdition), 21, 22),
     cex= 2,
     las= 1, 
     xlab= "PC1",
     ylab= "PC2")
legend("bottomright",
       pt.bg= adjustcolor(colorRampPalette(Cc)(nrow(pca)), 0.5),
       legend= pca$cdition,
       bty= "n",
       pch= ifelse(grepl("allograft", pca$cdition), 21, 22))
dev.off()

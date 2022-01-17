setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
dat <- readRDS("Rdata/RNA_tables_object.rds")

lapply(c("cutnrun", "allograft"), function(cdition){
    # PCA
    pdf(paste0("pdf/PCC_replicates/PCA_RNA_epiCancer_", cdition, ".pdf"),
        width = 7,
        height = 7)
    pca <- prcomp(as.matrix(dat$VSD_COUNTS[[cdition]], 1))
    pca <- as.data.table(pca$rotation, 
                         keep.rownames = T)
    pca[, cdition:= gsub("(.*)_.*", "\\1", rn)]
    plot(pca$PC1, 
         pca$PC2, col= factor(pca$cdition),
         pch= 19)
    legend("topright",
           legend= unique(pca$cdition),
           col= factor(unique(pca$cdition)),
           pch= 19,
           bty= "n")
    dev.off()
    
    # PCC
    pdf(paste0("pdf/PCC_replicates/RNA_epiCancer_", cdition, ".pdf"),
        width = 7,
        height = 7)
    vl_heatmap(cor(log2(as.matrix(dat$COUNTS[[cdition]], 1)+1)))
    dev.off()
    
    # # Dendrograms
    # dd <- as.dist(1 - cor(mat, method= "pearson"))
    # hc <- hclust(dd, method = "ward.D2")
    # pdf("pdf/PCC_replicates/RNA_epiCancer_Parreno_dendro.pdf", 
    #     width = 10, 
    #     height = 5)
    # par(bg= "grey")
    # plot(hc, 
    #      hang= -1, 
    #      cex= 0.6, 
    #      las= 1, 
    #      xlab= "Samples replicates",
    #      lend= 2)
    # points(seq(hc$order),
    #        rep(0, length(hc$order)),
    #        col= leg$Cc[hc$order],
    #        pch= leg$pch[hc$order],
    #        bg= "grey",
    #        cex= 0.6)
})



setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(scipen= 5)
require(data.table)
require(kohonen)
require(circlize)
require(pheatmap)

#------------------------------------------------------------#
# 1- Plot map
#------------------------------------------------------------#
som <- readRDS("Rdata/som_genes_FC_PH_vs_W.rds")
dat <- readRDS("Rdata/genes_features_final.rds")
dat <- na.omit(dat)
dat <- dat[, lapply(.SD, mean), som_cl, .SDcols= som_kcl:TSS_count_motif_cl99]
mdat <- melt(dat, id.vars = c("som_cl", "som_kcl"))
mdat <- mdat[order(som_cl)]
mdat[, cor1:= cor.test(value, som$codes[[1]][, "log2FoldChange_PHJ9_vs_WKD"])$estimate, variable]
mdat[, cor2:= cor.test(value, som$codes[[1]][, "log2FoldChange_PHJ11_vs_WKD"])$estimate, variable]
mdat[, cor3:= cor.test(value, som$codes[[1]][, "log2FoldChange_PH29_vs_W29"])$estimate, variable]
mdat[, cor:= max(abs(c(cor1, cor2, cor3))), variable]
mdat <- mdat[order(cor, -som_cl, decreasing= T)]
mdat <- mdat[!grepl("^padj", variable)]

size <- 28
pos <- matrix(seq(size*size), ncol= size, byrow = T)
pos <- pos[, c(7:size, 1:6)]
# pos <- pos[c(9:size, 1:8),]
# pos <- unlist(pos)
ord_kcl <- som$kcl[pos]

Cc1 <- colorRampPalette(c("cornflowerblue", "darkblue", "tomato", "lightblue", "limegreen", "gold", "white", "lightgrey"))

pdf("pdf/som_genes_FC_PH_vs_W_all_maps.pdf", 20, 20)
par(mfrow= c(5,5))
plot(som, "property", property= ord_kcl, palette= Cc1, shape= "straight", border= NA, main= "clusters")
add.cluster.boundaries(som, ord_kcl, lwd= 1, xpd= T)
# mdat[variable%in% unique(variable)[1:3],
mdat[,
    {
      lim <- boxplot(value, plot= F)$stats[c(1,5), 1]
      value[value<lim[1]] <- lim[1]
      value[value>lim[2]] <- lim[2]
      if(min(lim)<0)
      {
        Cc2 <- colorRamp2(c(lim[1], 0, lim[2]), c("cornflowerblue", "white", "tomato"))
      }else
      {
        Cc2 <- colorRamp2(lim, c("blue", "yellow"))
      }
      Cc2 <- colorRampPalette(Cc2(seq(lim[1], lim[2], length.out = 100)))
      plot(som, "property", property= value[pos], palette= Cc2, shape= "straight", zlim= lim, border= NA, main= variable)
      add.cluster.boundaries(som, ord_kcl, lwd= 1, xpd= T)
    }, variable]
dev.off()

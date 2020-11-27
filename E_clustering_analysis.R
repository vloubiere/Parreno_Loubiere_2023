setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(scipen= 5)
require(data.table)
require(org.Dm.eg.db)
require(kohonen)
require(circlize)
require(pheatmap)
require(clusterProfiler)

#------------------------------------------------------------#
# 1- Train Kohonen
#------------------------------------------------------------#
if(!file.exists("Rdata/som_genes_all.rds"))
{
  # Differential expression
  diff <- data.table(file= list.files("db/FC_tables_all_transcriptomes/", ".txt", full.names = T))
  diff[, exp := gsub("_FC.txt", "", basename(file))]
  setkeyv(diff, "exp")
  diff <- diff[c("PCXT1092A_vs_2A", "PSCSUZ21B842D_vs_42D", "EZ7312A_vs_2A", "SUZ1212A_vs_2A",
                 "PH18_vs_PH29", "PHJ9_vs_PH29", "PHJ11_vs_PH29",
                 "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18",
                 "PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29", "PHRNAI_vs_WRNAI",
                 "WTE1416_vs_120hED", "72hED_vs_120hED", "96hED_vs_120hED",
                 "72hED_vs_WTE1416", "96hED_vs_72hED", "120hED_vs_96hED")]
  diff <- diff[, fread(file), (diff)]
  # filter out genes that change at 18 degrees
  diff <- diff[!rn %in% diff[exp=="PH18_vs_W18" & padj<0.01, rn]]
  # Basemean
  baseMean <- diff[exp=="PH29_vs_W29", .(rn, baseMean_PH29_vs_W29= log2(baseMean+1))]
  diff <- dcast(diff, rn~exp, value.var = "log2FoldChange")
  # ChIP
  ChIP <- readRDS("Rdata/ChIP_normalized_enrichment_TSS_REs.rds")
  # Merge
  dat <- merge(diff, ChIP, by.x= "rn", by.y="FBgn")
  dat <- merge(dat, baseMean, by= "rn")
  dat <- melt(dat, id.vars = c("rn", "symbol"))
  # Clip outliers
  dat[!is.na(value), value:= 
      {
        lim <- quantile(value, c(0.01, 0.99), na.rm= T)
        value[value<lim[1]] <- lim[1]
        value[value>lim[2]] <- lim[2]
        value
      }, variable]
  dat <- dcast(dat, rn+symbol~variable, value.var = "value")
  
  # RNA
  RNAi <- dat[, .(PHJ9_vs_WKD, PHJ11_vs_WKD, PH29_vs_W29,
                  PH18_vs_PH29, PHJ9_vs_PH29, PHJ11_vs_PH29,
                  PHJ9_vs_PH18, PHJ11_vs_PH18, PH29_vs_PH18)]
  exp <- dat[, .(baseMean_PH29_vs_W29)]
  dev <- dat[, .(`72hED_vs_WTE1416`, `96hED_vs_72hED`, `120hED_vs_96hED`)]
  # mut <- dat[, .(PCXT1092A_vs_2A, PSCSUZ21B842D_vs_42D, SUZ1212A_vs_2A)]
  
  # ChIP
  TSS <- dat[, .(TSS_enr_PC, TSS_enr_PH, TSS_enr_PSC, TSS_enr_SUZ12, TSS_enr_POLIIGSE112868, 
                 TSS_enr_EYGSE112868, TSS_enr_H3K4me3, TSS_enr_H3K4me2, TSS_enr_H3K36me3, 
                 TSS_enr_H4K20me1, TSS_enr_H3K27ac, TSS_enr_H3K4me1, TSS_enr_H3K27me1, TSS_enr_H3K27me3, 
                 TSS_enr_H2AK118Ub, TSS_enr_H3K27me2)]
  RE <- dat[, .(RE_enr_PC, RE_enr_PH, RE_enr_PSC, RE_enr_SUZ12, RE_enr_POLIIGSE112868, 
                RE_enr_EYGSE112868, RE_enr_H3K4me3, RE_enr_H3K4me2, RE_enr_H3K36me3, 
                RE_enr_H4K20me1, RE_enr_H3K27ac, RE_enr_H3K4me1, RE_enr_H3K27me1, RE_enr_H3K27me3, 
                RE_enr_H2AK118Ub, RE_enr_H3K27me2)]
  
  # dat
  train <- list(RNAi, 
                do.call(cbind, list(exp, dev)),
                do.call(cbind, list(TSS, RE)))
  train <- lapply(train, function(x)
  {
    x <- as.matrix(x)
    rownames(x) <- dat[, paste0(rn, "__", symbol)]
    return(x)
  })
  
  mygrid <- somgrid(xdim= 28, ydim= 28, topo = 'hexagonal', toroidal = T)
  set.seed(1)
  som.model <- supersom(train, user.weights = c(1,5,5), grid = mygrid, rlen = 500)
  saveRDS(som.model, "Rdata/som_genes_TSS_RE.rds")
}

#------------------------------------------------------------#
# 2- Plot kohonen maps
#------------------------------------------------------------#
# Nodes clustering
som <- readRDS("Rdata/som_genes_TSS_RE.rds")
codes <- as.data.table(som$codes)
sub <- codes[, PHJ9_vs_WKD:PH29_vs_PH18]
set.seed(1)
sub[, kcl:= kmeans(.SD, 7)$cluster]
# sub[, kcl:= cutree(hclust(dist(.SD)), 8)]
# goi <- apply(sub, 1, function(x) any(abs(x)>1))
# sub[(goi), kcl:= kmeans(.SD, 6)$cluster]
# sub[is.na(kcl), kcl:= 7]

som_dat <- as.data.table(som$data)
som_dat[, cl:= som$unit.classif]
som_dat <- som_dat[, lapply(.SD, mean, na.rm= T), keyby= cl]
som_dat <- melt(som_dat, id.vars = "cl")
som_dat[grepl("_vs_", variable), c("zmin", "zmax"):= .(-3, 3)]
som_dat[grepl("baseMean_", variable), c("zmin", "zmax"):= .(0, 12)]
som_dat[grepl("_enr_", variable), c("zmin", "zmax"):= .(-2, 2)]

# Add genes groups
go <- data.table(file= list.files("db/GO_FBgn", full.names = T))
go[, variable:= gsub("_FlyBase_IDs.txt", "", basename(file))]
go <- go[, .(FBgn= unlist(tstrsplit(rownames(as.data.frame(som$data)), "__", keep=1)), 
             cl = som.model$unit.classif), (go)]
go[, value:= ifelse(FBgn %in% fread(file, header=F)$V1, 1, 0), file]
go <- go[, .(value= log2(sum(value)+1), zmin= 0, zmax= 2), keyby= .(variable, cl)] 
som_dat <- rbind(som_dat, go)

# PLOT kohonen
mat <- 1:54
j <- 55
for(i in c(2,3,4,5,9,10,14,15,19,20,25,29,30,35,38,39,40,55,58,59,60))
{
  mat <- c(mat[1:(i-1)], j, mat[i:length(mat)])
  j <- j+1
}

pdf("pdf/som_genes_TSS_RE.pdf", 20, 45)
layout(matrix(mat, byrow = T, ncol = 5))
Cc <- c("royalblue2","darkblue", "cornflowerblue","tomato","gold","lightgrey","lightblue")
Cc <- colorRampPalette(Cc)
plot(som, "property", property= sub$kcl, palette= Cc, shape= "straight", border= NA, main= "clusters")
som_dat[, 
        {
          value[value<zmin] <- zmin
          value[value>zmax] <- zmax
          Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))
          plot(som, "property", property= value, palette= Cc, shape= "straight", zlim= c(zmin, zmax), border= NA, main= variable)
          add.cluster.boundaries(som, sub$kcl, lwd= 1)
        }, .(variable, zmin, zmax)]
dev.off()

#------------------------------------------------------------#
# 3- Heatmap cluster
#------------------------------------------------------------#
# Heatmap clusters
hcl_dat <- as.data.table(som$data)
hcl_dat[, cl:= som$unit.classif]
hcl_dat[, kcl:= sub[cl, kcl]]
hcl_dat[, kcl_name:= paste("cluster", kcl, "(n=", .N, ")"), kcl]
h_sel <- c("PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29", "PH18_vs_PH29", "PHJ9_vs_PH29", 
         "PHJ11_vs_PH29", "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18")
h_dmat <- hcl_dat[, lapply(.SD, mean, na.rm= T), kcl_name, .SDcols= h_sel]

# Plot
Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))(100)
pheatmap(as.matrix(h_dmat, 1), col_labels = T, color= Cc, display_numbers = T, cluster_cols = F, 
         filename = "pdf/heatmap_som_clusters.pdf", height = 4, width = 6)

#------------------------------------------------------------#
# 4- GO
#------------------------------------------------------------#
if(!file.exists("Rdata/kcl_som_GO.rds"))
{
  go_dat <- as.data.table(som$data, keep.rownames = T)
  go_dat[, kcl:= sub$kcl[som$unit.classif]]
  go_dat[, c("fbgn", "symbol"):= tstrsplit(rn, "__")]
  entrez <- as.data.table(bitr(go_dat$fbgn, fromType= "FLYBASE", toType= "ENTREZID", OrgDb = org.Dm.eg.db))
  go_dat <- merge(go_dat[, .(fbgn, symbol, kcl)], entrez, by.x="fbgn", by.y="FLYBASE")
  
  res <- list()
  for(i in 1:6)
  {
    res[[i]] <- enrichGO(gene = go_dat[kcl==i, ENTREZID], OrgDb = org.Dm.eg.db, keyType = "ENTREZID", 
                           ont= "ALL", qvalueCutoff=0.01)
  }
  
  list.genes <- split(go_dat[kcl %in% 1:6, ENTREZID], go_dat[kcl %in% 1:6, kcl])
  res[[7]] <- compareCluster(list.genes, fun = "enrichGO", OrgDb= org.Dm.eg.db)
  saveRDS(res, "Rdata/kcl_som_GO.rds")
}
res <- readRDS("Rdata/kcl_som_GO.rds")

pdf("pdf/GO_kclusters_som.pdf", width = 14)
dotplot(res[[1]], showCategory= 20, title= "cluster 1")
dotplot(res[[2]], showCategory= 20, title= "cluster 2")
dotplot(res[[3]], showCategory= 20, title= "cluster 3")
dotplot(res[[4]], showCategory= 20, title= "cluster 4")
dotplot(res[[5]], showCategory= 20, title= "cluster 5")
dotplot(res[[6]], showCategory= 20, title= "cluster 6")
dotplot(res[[7]], showCategory= 7, title= "cluster 7")
dev.off()







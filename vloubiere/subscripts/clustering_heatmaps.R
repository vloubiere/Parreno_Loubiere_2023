# Clustering heatmaps
dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")

# Add PRC1 clusters
PRC1 <- get(load("external_data/SA2020_cl.list"))
PRC1 <- data.table(symbol= unlist(PRC1$genes_simp),
                   cl= rep(names(PRC1$genes_simp), lengths(PRC1$genes_simp)))
dat[PRC1, PRC1:= i.cl, on= "symbol"]
dat[is.na(PRC1), PRC1:= "UNBOUND"]
dat[, PRC1:= factor(PRC1, levels = c("aTSS A+B", "Enhancer A+B", "Polycomb A+B", "Repressed", "UNBOUND"))]
dat[, percPRC1:= .N/length(which(dat$cl==ccl))*100, .(ccl= cl, PRC1)]

# Add tissue-specificty score
tspe <- fread("external_data/tissue_specificty_score.txt")
dat[tspe, tissue_spe:= i.tissue_spe, on= "FBgn"]
dat[, tissue_spe_agg:= mean(unique(.SD[, .(FBgn, tissue_spe)])$tissue_spe, na.rm= T), cl]

# Add ChIP-Seq data
add <- readRDS("Rdata/available_data.rds")
add[dat, c("cl", "symbol"):= .(i.cl, i.symbol), on= "FBgn"]
add <- add[!is.na(cl)]

# Import motifs
mot <- readRDS("Rdata/selected_motifs.rds")

# Import paro
Paro <- data.table(file= list.files("db/FC_tables/", "Paro", full.names= T))
Paro <- Paro[, fread(file), file]
Paro[dat, cl:= i.cl, on= "V1==FBgn"]
Paro[, cdition:= gsub("RNA_Paro_2018_RNA_|_vs_RNA_WTED_FC.txt", "", basename(file))]

#---------------------------------------------------------------#

pdf("pdf/clustering_transcriptomes.pdf", width = 55, height = 7)
layout(matrix(1:12, ncol= 12), widths = c(0.7,0.7,0.3,0.7,0.7,0.3,0.4,0.4,1.5,3,3,0.3))
Cc <- c("blue", "cornflowerblue", "white", "tomato", "red")

par(xaxs= "i", 
    yaxs= "i", 
    las= 2,
    mar= c(12.1,8.1,8.1,6.1))

# Clustering heatmap #####
mat <- as.matrix(dcast(dat, 
                       cl+FBgn~cdition, 
                       value.var = "clipped_log2FoldChange")[,-2], 1)

vl_heatmap(mat, 
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "log2FC")
mtext("Transcriptome clustering", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)
box()
cllines <- 1-cumsum(dat[, .N, keyby= cl][,N])/nrow(dat)
abline(h= cllines[-length(cllines)])
axis(2, 
     labels = dat[, cl, keyby= cl][, cl],
     at = cllines-diff(c(1, cllines))/2, 
     tick = 0)

# PRC1_binding #####
.PRC1 <- dcast(dat, cl+FBgn~PRC1, fun.aggregate = function(x) ifelse(length(x)>0, 1, 0))
.PRC1 <- as.matrix(.PRC1[, `aTSS A+B`:UNBOUND])
vl_heatmap(mat, 
           newdata = .PRC1,
           col = c("white", "red"), 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "")
mtext("PRC1 binding (binary)", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)
box()
abline(h= cllines[-length(cllines)])
axis(2,
     labels = dat[, cl, keyby= cl][, cl],
     at = cllines-diff(c(1, cllines))/2,
     tick = 0)

# tissue specificity #####
par(mar= c(12.1,3.1,8.1,5.1))
.tspe <- as.matrix(dat[order(cl, FBgn), unique(tissue_spe), .(cl, FBgn)]$V1)
colnames(.tspe) <- "tissue_specifity_score"
vl_heatmap(mat, 
           newdata = .tspe,
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "")
mtext("tissue\nspecificty\nindex\n(high= hk, low= dev)", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)
box()
abline(h= cllines[-length(cllines)])
axis(2,
     labels = dat[, cl, keyby= cl][, cl],
     at = cllines-diff(c(1, cllines))/2,
     tick = 0)

# Plot transcriptome aggregates #####
par(mar= c(12.1,8.1,8.1,6.1))
agg <- as.matrix(dcast(dat, 
                       cl~cdition, 
                       fun.aggregate = mean,
                       value.var = "clipped_log2FoldChange"), 1)
vl_heatmap(agg, 
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "Mean\nlog2FC")
mtext("Transcriptome aggregate\n(/cluster)", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)
box()
axis(2,
     labels = dat[, paste0("cluster ", cl, " (", length(unique(FBgn)), ")"), keyby= cl][, V1],
     at = 1-dat[, cl, keyby= cl][, cl]/length(unique(dat$cl))+1/(2*length(unique(dat$cl))),
     tick = 0)

# PRC1_binding aggregate #####
.PRC1agg <- dcast(unique(dat[, .(cl, PRC1, percPRC1)]), cl~PRC1, value.var= "percPRC1")
.PRC1agg <- as.matrix(.PRC1agg[, -1])
.PRC1agg[is.na(.PRC1agg)] <- 0

vl_heatmap(agg, 
           newdata = .PRC1agg,
           col = c("white", "red", "darkred"), 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "Genes\npercentage", 
           display_numbers = T)
mtext("PRC1 binding\naggregate", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)
box()

# tissue speficity aggregate ######
.tspeagg <- as.matrix(unique(dat[, .(cl, tissue_spe_agg)]), 1)

par(mar= c(12.1,3.1,8.1,5.1))
vl_heatmap(agg, 
           newdata = .tspeagg,
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           show_rownames = F, 
           legend_title = "", 
           display_numbers = T)
mtext("tissue\nspecificty\nindex\naggregate\n(high= hk, low= dev)", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)
box()

# Plot dev transcriptomes ######
transcriptomes <- as.matrix(dcast(add[.id=="Transcriptome"], 
                                  cl~cdition, 
                                  fun.aggregate = mean, 
                                  na.rm=T,
                                  value.var = "score"), 1)

par(mar= c(12.1,2.1,8.1,6.1))
vl_heatmap(mat = agg, 
           newdata = transcriptomes[, c(2,3,1)],
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= F,
           breaks= c(-3, 0, 3),
           legend_pos = c(0.8, 1.15, 1, 1.27),
           legend_title = "Mean\nlog2FC")
mtext("ED dev.\ntranscriptomes", 
      side= 3, 
      cex= 0.5, 
      las=1, 
      line = 1)

# Paro ecdysone transcriptomes ######
.paro <- as.matrix(dcast(Paro[!is.na(cl)], cl~cdition, value.var = "log2FoldChange", fun.aggregate = mean), 1)
vl_heatmap(.paro, 
           cluster_rows = F, 
           display_numbers = T, 
           breaks = c(-8,0,8))

# Plot ChIP aggregate ######
ChIP <- as.matrix(dcast(add[.id=="ChIP"], 
                        cl~cdition, 
                        fun.aggregate = mean, 
                        na.rm=T,
                        value.var = "score"), 1)
colnames(ChIP)[colnames(ChIP)=="RPKM"] <- "RPKM wt ED" 

par(mar= c(12.1,2.1,8.1,4.1))
vl_heatmap(mat = agg, 
           newdata = scale(ChIP),
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= T, 
           breaks = c(-4,0,4),
           legend_pos = c(0.8, 1.02, 1, 1.035),
           legend_title = "z score")
mtext("ChIP-Seq enrichment aggregate", 
      side= 3, 
      cex= 0.7, 
      las=1, 
      line = 4)

# Plot motif enrichment ALL MOTIFS #####
motifs <- as.matrix(dcast(mot$all_REs,
                          cl~symbol,
                          na.rm=T,
                          value.var = "log2FE"), 1)
vl_heatmap(mat = agg, 
           newdata = motifs,
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= T,
           breaks= c(-3, 0, 3),
           legend_pos = c(0.8, 1.017, 1, 1.025),
           legend_title = "odd\nratio")
mtext("TF motifs enrichment all REs", 
      side= 3, 
      cex= 0.7, 
      las=1, 
      line = 4)

# Plot motif enrichment TSSs ######
motifs <- as.matrix(dcast(mot$TSSs,
                          cl~symbol,
                          na.rm=T,
                          value.var = "log2FE"), 1)
vl_heatmap(mat = agg, 
           newdata = motifs,
           col = Cc, 
           cluster_rows= F, 
           cluster_cols= T,
           breaks= c(-3, 0, 3),
           legend_pos = c(0.8, 1.017, 1, 1.025),
           legend_title = "odd\nratio")
mtext("TF motifs enrichment TSSs", 
      side= 3, 
      cex= 0.7, 
      las=1, 
      line = 4)

# Plot motif enrichment Enhancers ######
motifs <- as.matrix(dcast(mot$enhancers,
                          cl~symbol,
                          na.rm=T,
                          value.var = "log2FE"), 1)

par(mar= c(12.1,3.1,8.1,5.1))
vl_heatmap(mat = agg, 
           newdata = motifs,
           col = Cc, 
           cluster_rows= F,
           breaks= c(-3, 0, 3),
           legend_pos = c(0.8, 1.2, 1, 1.3),
           legend_title = "odd\nratio")

mtext("TF motifs \nenrichment REs", 
      side= 3, 
      cex= 0.7, 
      las=1, 
      line = 4)

dev.off()


setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
source("/groups/stark/vloubiere/scripts/R_functions/bedtools_wrappers_1.0.R")
require(circlize)
require(pheatmap)
require(data.table)
require(gridExtra)
require(clusterProfiler)
require(org.Dm.eg.db)
require(rtracklayer)
require(colorRamps)
require(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
require(motifmatchr)
require(PWMEnrich) # module load gsl/2.1-foss-2017a
require(TFBSTools)
require(seqLogo)

#--------------------------------------------------------------------#
# 1- Import and annotate data
#--------------------------------------------------------------------#
dat <- data.table(file= list.files("FC_tables", "PH18_vs_W18|PH29_vs_W29|PHJ9_vs_WKD|PHJ11_vs_WKD", full.names = T))
dat[, exp:= gsub("_FC.txt", "", basename(file))]
dat <- dat[, fread(file), c(colnames(dat))]

# Adcc chromatin types (overlapping the largest fraction of the gene body)
tss <- genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
tss$name <-  mapIds(org.Dm.eg.db, key= tss$gene_id, column="SYMBOL", keytype="FLYBASE", multiVals="first")
chrom <- import("chrom_types/chromatin_types_SA_2020_dm3.bed")
t1 <- tempfile(fileext = ".bed")
t2 <- tempfile(fileext = ".bed")
export(tss, t1)
export(chrom, t2)
chrmm <- fread(cmd= paste("/software/2020/software/bedtools/2.27.1-foss-2018b/bin/intersectBed -wo -a", t1, "-b", t2))
chrmm <- chrmm[, .(V14= sum(V13)), .(V4, V10)]
chrmm <- chrmm[, V10[which.max(V14)], V4]
dat[chrmm, chrmm:= i.V1, on= "rn==V4"]

# PH binding
PH_target <- get(load("/groups/stark/vloubiere/projects/SA_2019/data_SA_2019/objects/cl.list"))$genes
PH_target <- data.table(symbol= unlist(PH_target), type= rep(names(PH_target), lengths(PH_target)))
PH_target <- PH_target[!is.na(type)]
PH_target[, N:= .N, .(symbol, type)]
PH_target <- PH_target[, .(type= type[which.max(N)]), symbol]
dat[PH_target, PH_target:= i.type, on= "rn==symbol"]

#--------------------------------------------------------------------#
# 2- wRNAi vs phRNAi escapers G.O
#--------------------------------------------------------------------#
c_genes <- dat[exp=="PH18_vs_W18" & padj<0.01]

pdf("pdf/list_wRNAi_phRNAi_18_escaper_genes.pdf", width = 7, 18)
plot.new()
grid.table(c_genes[log2FoldChange>0, .(rn, log2FoldChange)][order(log2FoldChange, decreasing = T)])
mtext("Up-regulated phRNAi vs wRNAi 18C")
plot.new()
grid.table(c_genes[log2FoldChange<0, .(rn, log2FoldChange)][order(log2FoldChange)])
mtext("Down-regulated phRNAi vs wRNAi 18C")
dev.off()

# # Gene ontologies
# go_up <- enrichGO(bitr(c_genes[log2FoldChange>0, rn], fromType= "SYMBOL", toType= "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID, OrgDb = org.Dm.eg.db)
# dotplot(go_up, showCategory= 10)
# go_down <- enrichGO(bitr(c_genes[log2FoldChange<0, rn], fromType= "SYMBOL", toType= "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID, OrgDb = org.Dm.eg.db)
# dotplot(go_down, showCategory= 10)
# go_all <- enrichGO(bitr(c_genes$rn, fromType= "SYMBOL", toType= "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID, OrgDb = org.Dm.eg.db)
# dotplot(go_all, showCategory= 10)

#--------------------------------------------------------------------#
# 3- Up and down heatmap
#--------------------------------------------------------------------#
sub <- copy(dat)
sub[padj<0.01 & abs(log2FoldChange)>1, binFC:= ifelse(log2FoldChange>0, 1, -1)]
sub <- sub[is.na(binFC), binFC:=0][binFC!=0]
dmat <- dcast(sub, rn+chrmm+PH_target~exp, value.var = "binFC", fill = 0)
# Clustering based on changes
dmat[PH18_vs_W18!=0, cl:= 1]
dmat[PH18_vs_W18==0, cl:= .GRP+1, .(PHJ9_vs_WKD, PHJ11_vs_WKD, PH29_vs_W29)]
dmat[, N:= .N, cl]
dmat[order(N, decreasing= T), cl:= .GRP, cl]
dmat <- dmat[order(cl)]
# format matrix
mat <- as.matrix(dmat[, c(1,4,7,6,5)], 1)
# Annot
annot <- data.frame(cluster= formatC(dmat$cl, flag = "0", width= 2),
                    PH_target= dmat$PH_target,
                    chrom_type= dmat$chrmm,
                    row.names = dmat$rn)
# Annot colors
cl_col <- c(matlab.like2(11), grey.colors(10))
names(cl_col) <- unique(annot$cluster)
PH_target <- c("tomato", "brown3", "limegreen", "darkgreen", "#54BFED", "royalblue1", "#D3D3D3")
names(PH_target) <- c("aTSS A","aTSS B","Enhancer A","Enhancer B","Polycomb A","Polycomb B","Repressed")
chrom_type <- c("tomato", "brown3", "#74C27A", "#54BFED", "#D3D3D3")
names(chrom_type) <- c("aTSS", "aTTS", "Enhancer", "PcG", "Null")
annot_col <- list(cluster= cl_col, PH_target= PH_target, chrom_type= chrom_type)

# Heatmap
Cc <- c("cornflowerblue", "lightgrey", "tomato")
lb <- c(-1,0,1)
ll <- c("down-regulated", "stable", "up-regulated")
gr <- cumsum(table(dmat$cl))
pheatmap(mat, cluster_cols = F, cluster_rows = F, show_rownames = F, col= Cc, legend_breaks = lb, legend_labels = ll, 
         gaps_row = gr, annotation_row = annot, annotation_colors = annot_col, 
         filename = "pdf/heatmap_clustering_FC_patterns.pdf", width = 6, height = 8,
         cellwidth = 20, cellheight = 0.15)

dmat <- dmat[, .(symbol= rn, chromatin_type= chrmm, PH_target, cluster= cl)]
dmat[cluster==1, metacluster:= 1]
dmat[cluster %in% c(2,3,5), metacluster:= 2]
dmat[cluster==4, metacluster:= 3]
dmat[cluster %in% c(9, 11), metacluster:= 4]
dmat[cluster %in% c(6, 7), metacluster:= 5]
fwrite(dmat, "data/gene_clusters.txt", col.names = T, row.names = F, sep= "\t", quote= F, na = NA)

#--------------------------------------------------------------------#
# 4- metaclusters analysis
#--------------------------------------------------------------------#
sub <- fread("data/gene_clusters.txt")
setkeyv(dat, "rn")

pdf("pdf/metaclusters_profiles.pdf", width = 11, height = 8)
par(mfrow= c(2, 3), mar= c(8,5,5,5))
for(i in 1:5)
{
  current <- dcast(dat[sub[metacluster==i, symbol]], rn~exp, value.var = "log2FoldChange")
  plot(NA, xlim= c(1,4), ylim= c(-6, 6), ylab= "log2FoldChange", xaxt= "n", xlab= "", las = 2)
  current[, lines(1:4, c(PH18_vs_W18, PHJ9_vs_WKD, PHJ11_vs_WKD, PH29_vs_W29), col= adjustcolor("lightgrey", 0.6)), c(colnames(current))]
  abline(h= 0, lty= 2)
  axis(1, at= 1:4, c("PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29"), las= 2)
  mtext(paste0("metaCluster ", i), line= 0.5)
  mtext(c("Rescuable down", "Stably down", "Rescuable up", "Stable up", "Transciently up")[i], line= 2)
}
dev.off()

#--------------------------------------------------------------------#
# 5- chromatin types
#--------------------------------------------------------------------#
sub <- fread("data/gene_clusters.txt")
perc <- sub[!is.na(metacluster), .(perc= as.numeric(.N)), .(chromatin_type, metacluster)]
perc[, perc:= perc/sum(perc)*100, metacluster]
perc <- rbind(perc, dat[, .(perc= .N/nrow(dat)*100, metacluster= "genome"), .(chromatin_type= chrmm)])
res <- as.matrix(dcast(perc, chromatin_type~metacluster, value.var = "perc"), 1)[c(5,6,2,4,3,1),]
Cc <- c("tomato", "brown3", "#74C27A", "#54BFED", "#D3D3D3", "white")

pdf("pdf/barplot_chromatin_types.pdf")
barplot(res, las= 1, col= Cc, ylim= c(0, 150), ylab= "Percentage of genes")
legend("topright", fill = Cc, legend = c("aTSS", "aTTS", "Enhancer", "PcG", "Null", "NA"), bty= "n")
dev.off()

#--------------------------------------------------------------------#
# 6- Gene ontologies
#--------------------------------------------------------------------#
# GO comparecluster
sub <- fread("data/gene_clusters.txt")

# gene.list <- split(sub$symbol, sub$metacluster)
# gene.list <- lapply(gene.list, function(x) bitr(x, fromType= "SYMBOL", toType= "ENTREZID", OrgDb = org.Dm.eg.db)$ENTREZID)
# res <- compareCluster(gene.list, fun = "enrichGO", OrgDb= org.Dm.eg.db)
# # Plot
# pdf(paste0("pdf/GO_enrichment_metaclusters.pdf"), 20, 12)
# dotplot(res, showCategory= 10, title= "GO enrich metaclusters")
# dev.off()

# single <- enrichGO(gene.list[[5]], OrgDb = org.Dm.eg.db)
# dotplot(single, showCategory= 10, title= "GO enrich metaclusters")

#--------------------------------------------------------------------#
# 7- Motifs at promoters
#--------------------------------------------------------------------#
# ### Motifs counts  ############################
# # load motifs and select the ones with Dmel names
# load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
# DT <- as.data.table(TF_clusters_PWMs$metadata)
# sel <- DT[!is.na(Dmel), motif_name]
# sel <- match(sel, name(TF_clusters_PWMs$All_pwms_log_odds))
# 
# # Count hits at TSSs
# tss <- resize(genes(TxDb.Dmelanogaster.UCSC.dm3.ensGene), 1, "start")
# tss <- resize(tss, 250, "center")
# tss$name <-  mapIds(org.Dm.eg.db, key= tss$gene_id, column="SYMBOL", keytype="FLYBASE", multiVals="first")
# hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], tss, genome= "dm3", p.cutoff= 1e-4, bg="even", out= "scores")
# 
# counts <- as.matrix(motifCounts(hit))
# colnames(counts) <- DT[name(TF_clusters_PWMs$All_pwms_log_odds[sel]), paste0(Dmel, "_motif_count.", motif_name), on= "motif_name"]
# rownames(counts) <- tss$name
# counts <- as.data.table(counts, keep.rownames = T)
# 
# # Add metaclusters
# counts[sub, metacluster:= i.metacluster, on= "rn==symbol"]
# counts <- melt(counts, id.vars = c("rn", "metacluster"))
# counts[is.na(metacluster), metacluster:= 0]
# 
# # Collapse on prot name based on counts
# counts[, Dmel:= tstrsplit(variable, "_motif_count", keep= 1), variable]
# counts[, check := sum(value), .(Dmel, variable)]
# counts[, check2:= check==max(check), .(Dmel)]
# counts <- counts[(check2), !c("check", "check2")]
# 
# # Compute fisher tests
# counts[, c("fisher_OR_1", "fisher_pval_1"):= fisher.test(table(metacluster==1, value>0))[c("estimate", "p.value")], Dmel]
# counts[, c("fisher_OR_2", "fisher_pval_2"):= fisher.test(table(metacluster==2, value>0))[c("estimate", "p.value")], Dmel]
# counts[, c("fisher_OR_3", "fisher_pval_3"):= fisher.test(table(metacluster==3, value>0))[c("estimate", "p.value")], Dmel]
# counts[, c("fisher_OR_4", "fisher_pval_4"):= fisher.test(table(metacluster==4, value>0))[c("estimate", "p.value")], Dmel]
# counts[, c("fisher_OR_5", "fisher_pval_5"):= fisher.test(table(metacluster==5, value>0))[c("estimate", "p.value")], Dmel]
# cols <- c("Dmel", grep("^fisher", colnames(counts), value= T))
# res <- unique(counts[, ..cols])
# saveRDS(res, "data/motif_counts_collapsed_1.0.rds")

counts <- readRDS("data/motif_counts_collapsed_1.0.rds")
cols <- gsub("pval", "padj", grep("pval", colnames(counts), value= T))
counts[, (cols):= lapply(.SD, p.adjust), Dmel, .SDcols= patterns("pval")]
counts <- counts[fisher_padj_1<0.01|fisher_padj_2<0.01|fisher_padj_3<0.01|fisher_padj_4<0.01|fisher_padj_5<0.01]
counts <- counts[fisher_OR_1>2|fisher_OR_2>2|fisher_OR_3>2|fisher_OR_4>2|fisher_OR_5>2]
counts <- as.matrix(counts[, .(Dmel, fisher_OR_1, fisher_OR_2, fisher_OR_3, fisher_OR_4, fisher_OR_5)], 1)
colnames(counts) <- paste("metacl", 1:5)
rownames(counts) <- unlist(tstrsplit(rownames(counts), "_motif_count", keep= 1))
counts <- log2(counts)
counts[counts==Inf|counts==-Inf] <- NA
pheatmap(t(counts), cluster_rows = F, breaks= seq(-3, 3, length.out = 101), 
         col= colorRampPalette(c("cornflowerblue", "white", "tomato"))(100),
         filename = "pdf/motif_enrichment_TSSs.pdf", width = 15, cutree_cols = 7)

# # Circular version
# pdf("pdf/circular_heatmap.pdf")
# 
# circos.par(gap.after= c(rep(2, 20), 15)) # Set gap afeter last cluster (21 clusters)
# # Plot FC
# col_fun1 <- colorRamp2(c(-1, 0, 1), c("cornflowerblue", "lightgrey", "tomato"))
# circos.heatmap(mat, col = col_fun1, split = as.character(unlist(tstrsplit(rownames(mat), "__", keep=2))), track.height= 0.4)
# # Add colnames
# circos.track(track.index = 1, panel.fun = function(x, y) {
#   if(CELL_META$sector.numeric.index == 21) { # the last sector
#     cn = colnames(mat)
#     n = length(cn)
#     circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
#                 1:n - 0.5, cn, 
#                 cex = 0.5, adj = c(0, 0.5), facing = "inside")
#   }
# }, bg.border = NA)
# # Plot cluster names
# Cc <- matlab.like2(21)
# names(Cc) <- 1:21
# circos.heatmap(as.matrix(as.numeric(annot$cluster)), col = Cc, track.height= 0.1)
# circos.clear()
# dev.off()


# clmm_names <- c("aTSS","Enhancer","Polycomb","Intermediate")
# clmm_col <- c("tomato", "#74C27A", "#54BFED", "#D3D3D3")
# cl2_names <- c("aTSS A+B","Enhancer A+B","Polycomb A+B","Repressed")
# cl2_col <- c("#F27676", "#74C27A", "#54BFED", "#D3D3D3")
# cl_names <- c("aTSS A","aTSS B","Enhancer A","Enhancer B","Polycomb A","Polycomb B","Repressed")
# som_col <- c("tomato", "brown3", "#74C27A", "#54BFED", "#D3D3D3")
# som_names <- c("aTSS", "aTES", "Enhancer", "PcG", "Repressed")
# cl_col <- c("tomato", "brown3", "limegreen", "darkgreen", "#54BFED", "royalblue1", "#D3D3D3")






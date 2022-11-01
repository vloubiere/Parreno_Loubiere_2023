setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]
dat <- vl_resizeBed(dat, "start", 750, 250)
sel <- vl_Dmel_motifs_DB_full[Dmel %in% c("Stat92E", "kay/Jra")
                              | grepl("ZEB", motif_ID, ignore.case = T), .(motif_ID, motif_cluster, Dmel)]
sel[is.na(Dmel), Dmel:= "ZEB1"]
counts <- vl_motif_counts(dat[, .(seqnames, start, end)], 
                          sel= sel$motif_ID,
                          genome = "dm6")
dat <- cbind(dat[, .(recovery, symbol)], counts)
setnames(dat, 
         as.character(sel$motif_ID),
         make.unique(sel$Dmel))

pdf("pdf/recovery_promoters_motif_counts_heatmap.pdf", 5, 12)
par(mar= c(4,10,1,7),
    las= 2,
    cex.axis= 0.5)
vl_heatmap(as.matrix(dat[, symbol:ZEB1.2], 1), 
           cluster_rows= F,
           cluster_cols= F,
           row_clusters= dat$recovery,
           breaks= c(0,4), 
           col= c("blue", "yellow"), 
           legend_title = "Motif counts")
dev.off()

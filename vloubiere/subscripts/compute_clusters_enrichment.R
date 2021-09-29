require(TFBSTools)
require(motifmatchr)
require(vlfunctions)
dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")
feat <- readRDS("Rdata/ED_REs_features_final.rds")
feat[dat, cl:= i.cl, on="FBgn"]
feat[is.na(cl), cl:= 0]

# Compute mean FC
cl <- dcast(dat, 
            cl~cdition, 
            value.var= "log2FoldChange", 
            fun.aggregate = mean)

# N REs/gene
cl[feat[, .SD[, .N, FBgn], cl][, mean(N), cl], mean_N_REs:= i.V1, on= "cl"]

# Compute PRC1 bound fraction
cl[feat[, length(which(PRC1_bound))/.N, cl], PRC1_bound_fraction:= i.V1, on= "cl"]

# Compute mean CHIP
ChIP <- grep("_score$", colnames(feat), value = T)
ChIP <- feat[, lapply(.SD, median), cl, .SDcols= ChIP]
cl <- merge(cl, ChIP)

# Compute Enriched motifs
if(!file.exists("Rdata/fisher_tests_transcriptomes_clusters.rds"))
{
  mot <- grep("^motif:", colnames(feat), value = T)
  mot <- melt(feat, 
              id.vars= c("cl", "FBgn"),
              measure.vars = mot)
  mot <- mot[, .(value= sum(value)), .(FBgn, variable, cl)]
  mot <- mot[, {
    .var1 <- value>0
    .var2 <- cl
    .SD[cl!=0, {
      fisher.test(.var1, .var2==cl)[c("estimate", "p.value")]
    }, cl]
  }, variable]
  mot[, log2OR:= log2(estimate)]
  mot[, padj:= p.adjust(p.value, "fdr")]
  saveRDS(mot, 
          "Rdata/fisher_tests_transcriptomes_clusters.rds")
} 
mot <- readRDS("Rdata/fisher_tests_transcriptomes_clusters.rds")
sel <- mot[order(-abs(log2OR))][, .GRP, variable][GRP<=100, variable]
mot <- mot[variable %in% sel]
mot[, Dmel:= tstrsplit(variable, ":", keep= 3)]
mot <- mot[variable %in% mot[, variable[which.max(abs(log2OR))], Dmel]$V1]
mot <- dcast(mot, cl~Dmel, value.var = "log2OR")

mat <- as.matrix(mot, 1)
mat[mat==(-Inf)] <- min(mat[is.finite(mat)])
mat[mat==(Inf)] <- max(mat[is.finite(mat)])

vl_heatmap(mat, 
           cluster_rows = F, 
           breaks = c(-3, 0, 3))




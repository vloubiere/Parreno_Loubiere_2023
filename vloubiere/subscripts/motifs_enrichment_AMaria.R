# PRC1 peaks
PRC1 <- get(load("external_data/SA2020_cl.list"))
PRC1 <- cbind(as.data.table(PRC1$ED_summits), 
              data.table(cluster= PRC1$clusters_simp))
colnames(PRC1) <- c("seqnames", "start", "end", "ID", "pval", "cluster")

# ATAC-Seq control peaks
ctls <- get(load("external_data/SA2020_ctl.list"))
ctls <- cbind(as.data.table(ctls$ED_summits), 
              data.table(cluster= ctls$clusters_simp))
colnames(ctls) <- c("seqnames", "start", "end", "ID", "pval", "cluster")

# All regions
dat <- rbindlist(list("+"= PRC1, "-"= ctls), idcol = "PRC1_binding")
dat[, start:= start-124]
dat[, end:= end+125]

# count motifs
if(!file.exists("Rdata/counts_motifs_SA2020_clusters.rds"))
{
  hit <- matchMotifs(vl_Dmel_motifs_DB$All_pwms_log_odds, 
                     GRanges(dat[, seqnames:end], mcols= dat[, .(ID, cluster, PRC1_binding)]), 
                     genome= "dm6", 
                     p.cutoff= 5e-4, 
                     bg= "even", 
                     out= "scores")
  counts <- as.matrix(motifCounts(hit))
  colnames(counts) <- name(vl_Dmel_motifs_DB$All_pwms_log_odds)
  counts <- as.data.table(counts)
  res <- cbind(dat, counts)
  saveRDS(res, 
          "Rdata/counts_motifs_SA2020_clusters.rds")
}else
{
  res <- readRDS("Rdata/counts_motifs_SA2020_clusters.rds")
}

# Compute enrichment
cols <- colnames(res)[-c(1:7)]
.m <- melt(res, id.vars = c("PRC1_binding", "seqnames", "start", "end", "ID", "pval", "cluster"))
ftest <- .m[,
            {
              .c <- fisher.test(table(PRC1_binding=="+", value>0))
              .(pval= .c$p.value, est= .c$estimate)
            }, .(cluster, variable)]
ftest[, padj:= p.adjust(pval, method= "fdr"), cluster]
ftest[, odd_ratio:= log2(est)]
ftest[, check:= any(padj<0.00001) 
              & any(abs(odd_ratio)>1)
              & grepl("^bergman|^flyfactorsurvey|^homer|^jaspar|^idmmpmm", variable), variable]

mat <- dcast(ftest[(check)], variable~cluster, value.var = "odd_ratio")

pdf("pdf/AMaria_motifs_PRC1_clusters.pdf", height = 25, width = 7)
par(mar= c(3,25,3,5))
vl_heatmap(as.matrix(mat, 1), breaks = c(-2,0,2))
dev.off()

setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(glmnet)

# Import peaks ----
dat <- readRDS("Rdata/final_ATAC_table.rds")
dat <- dat[!is.na(log2FoldChange_PH29) & !is.na(log2FoldChange_PHD11)]
dat <- vl_resizeBed(dat,
                    "center",
                    250,
                    250,
                    genome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)

# Motif counts ----
if(!file.exists("db/counts/motifs/review_ATAC_peaks_motif_counts_500.rds"))
{
  # Select motifs whose cognate TF is transcribed
  sel <- vl_motifs_DB_v2[Dmel %in% readRDS("Rdata/final_gene_features_table.rds")$symbol, motif_ID]
  # Counts
  counts <- vl_motif_counts(dat,
                            genome= BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6,
                            sel= sel,
                            motifDB= vl_motifs_DB_v2)
  saveRDS(counts,
          "db/counts/motifs/review_ATAC_peaks_motif_counts_500.rds")
} else
  counts <- readRDS("db/counts/motifs/review_ATAC_peaks_motif_counts_500.rds")

# LASSO prediction ----
if(!file.exists("db/model/ATAC_FC_PH29_motif_LASSO.rds"))
{
  m29 <- vl_trainLASSO(response = dat$log2FoldChange_PH29,
                       counts= as.matrix(counts))
  saveRDS(m29,
          "db/model/ATAC_FC_PH29_motif_LASSO.rds")
}else
  m29 <- readRDS("db/model/ATAC_FC_PH29_motif_LASSO.rds")

if(!file.exists("db/model/ATAC_FC_PHD11_motif_LASSO.rds"))
{
  md11 <- vl_trainLASSO(response = dat$log2FoldChange_PHD11,
                        counts= as.matrix(counts))
  saveRDS(md11,
          "db/model/ATAC_FC_PHD11_motif_LASSO.rds")
}else
  md11 <- readRDS("db/model/ATAC_FC_PHD11_motif_LASSO.rds")

# Select relevant motifs ----
PH29 <- as.data.table(as.matrix(m29$beta), keep.rownames = T)
PH29[vl_motifs_DB_v2, name:= i.Dmel, on= "rn==motif_ID"]
PH29 <- PH29[, .SD[which.max(abs(s0))], name]
PH29 <- PH29[order(abs(s0), decreasing = T)]

PHD11  <- as.data.table(as.matrix(md11$beta), keep.rownames = T)
PHD11[vl_motifs_DB_v2, name:= i.Dmel, on= "rn==motif_ID"]
PHD11 <- PHD11[, .SD[which.max(abs(s0))], name]
PHD11 <- PHD11[order(abs(s0), decreasing = T)]

sel29 <- PH29[, .SD[which.max(abs(s0))], name]
seld11 <- PHD11[, .SD[which.max(abs(s0))], name]
sel <- rbind(sel29[1:25],
             PHD11[1:25])
sel <- sel[, .SD[which.max(abs(s0))], name]

# Select relevant motif counts and train lm ----
subCounts <- counts[, sel$rn, with= F]
setnames(subCounts, sel$name)
subCounts <- as.matrix(subCounts)
saveRDS(subCounts,
        "db/counts/motifs/review_selected_ATAC_peaks_motif_counts_500.rds")

# Linear models ----
PH29 <- lm(dat$log2FoldChange_PH29~subCounts)
PHD11 <- lm(dat$log2FoldChange_PHD11~subCounts)
sum29 <- as.data.table(summary(PH29)$coeff, keep.rownames = T)
sum29 <- sum29[-1]
setnames(sum29, names(sum29)[-1], function(x) paste0(x, ".29"))
sumd11 <- as.data.table(summary(PHD11)$coeff, keep.rownames = T)
sumd11 <- sumd11[-1]
setnames(sumd11, names(sumd11)[-1], function(x) paste0(x, ".d11"))
sum <- merge(sumd11, sum29)
sum[, rn:= gsub("subCounts", "", rn)]
setnames(sum, "rn", "name")
saveRDS(sum,
        "db/model/ATAC_FC_PH29_PHD11_motif_lm.rds")

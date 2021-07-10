
# I could not figure out how the qvalues are computed by novogene!

# dat <- readRDS("Rdata/DNA_novogene_full.rds")
# dat <- readRDS("Rdata/DNA_novogene_clean_final.rds")
# 
# # Example compute difference ph29_vs_ph29_t5
# diff <- dat[variable=="ph29_t5", .(CHROM, POS, ALT_UNIQ_ID, allele, ph29_t5_allele_counts= allele_counts)]
# diff[dat[variable=="ph29"], ph29_allele_counts:= i.allele_counts, on= "ALT_UNIQ_ID"]
# 
# dat[, allele_1_counts:= tstrsplit(allele_counts, ",", keep = idx), .(idx= as.numeric(gsub("(.*)/.*", "\\1", allele))+1)]
# dat[, allele_2_counts:= tstrsplit(allele_counts, ",", keep = idx), .(idx= as.numeric(gsub(".*/(.*)", "\\1", allele))+1)]
# dat[, allele_1_counts:= as.numeric(allele_1_counts)]
# dat[, allele_2_counts:= as.numeric(allele_2_counts)]
# 
# ph29 <- merge(dat[variable=="ph29", .(CHROM, POS, ALT_UNIQ_ID, allele_1_counts, allele_2_counts)],
#               dat[variable=="ph29_t5", .(CHROM, POS, ALT_UNIQ_ID, allele_1_counts, allele_2_counts)],
#               by= "ALT_UNIQ_ID")
# ph29[, c("pval", "log2OR"):= {
#   mat <- matrix(c(allele_1_counts.x, allele_2_counts.x, allele_1_counts.y, allele_2_counts.y), byrow = T, ncol = 2)
#   res <- fisher.test(mat)
#   .(res$p.value, log2(res$estimate))
# }, .(allele_1_counts.x, allele_2_counts.x, allele_1_counts.y, allele_2_counts.y)]
# ph29[, padj:= p.adjust(pval, method = "fdr")]

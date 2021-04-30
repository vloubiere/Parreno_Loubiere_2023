dat <- fread("../../genomes/dm6/gene_rpkm_report_fb_2021_02.tsv")
colnames(dat) <- gsub("#", "", colnames(dat))

clean <- dcast(dat[Parent_library_name == "modENCODE_mRNA-Seq_development"], FBgn+GeneSymbol~RNASource_name, value.var = "RPKM_value")
clean <- na.omit(clean[rowSums(clean[, 3:ncol(clean)])>1])
clean <- melt(clean, id.vars = c("FBgn", "GeneSymbol"), value.name = "RPKM", variable.name = "condition")
class <- unique(dat[, .(Parent_library_name, RNASource_name)])
clean[class, RNASeq:= i.Parent_library_name, on= "condition==RNASource_name"]

clean[, quant:= quantile(RPKM, 0.75, na.rm = T), condition]
clean[, class := ifelse(all(RPKM>quant), "housekeeping", "developmental"), FBgn]
clean[, tissue_spe := length(which(RPKM>quant))/length(RPKM), FBgn]

res <- unique(clean[, .(FBgn, GeneSymbol, class, tissue_spe)])
fwrite(res, "Rdata/tissue_specificty_score.txt", sep= "\t", quote=F, col.names = T, row.names = F)

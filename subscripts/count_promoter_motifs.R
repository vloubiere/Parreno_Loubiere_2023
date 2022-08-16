setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

cl <- rbindlist(list(cutnrun= cutnrun$regions,
                     allograft= allograft$regions), idcol = T)
FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
FBGN <- as.data.table(FBGN)
FBGN <- FBGN[type=="gene" & seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")]
promoters <- vl_resizeBed(FBGN, center= "start", upstream = 1000, downstream = 500)
promoters <- promoters[start>0, .(seqnames, start, end, strand, gene_id, gene_symbol)]
counts <- vl_motif_counts(bed = promoters[, seqnames:end], 
                          genome = "dm6")
rownames(counts) <- promoters[, paste0(seqnames, ":", start, "-", end, ":", strand, ":", gene_id, ":", gene_symbol)]
colnames(counts) <- vl_Dmel_motifs_DB_full[match(colnames(counts), vl_Dmel_motifs_DB_full$motif), uniqName_noSpecialChar]
saveRDS(counts, 
        "Rdata/motif_counts_promoters.rds")

setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(vlfunctions)
require(GenomicRanges)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat[, c("seqnames", "start"):= tstrsplit(id, ":|-", keep= c(1,2))]
dat[, start:= as.numeric(start)]
dat[, total:= .N, .(cdition, class, alt_class)]

bins <- vl_binBSgenome("dm6", 
                       bins_width = 1e5, 
                       restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4","chrX","chrY"))
# counts hits per bin
counts <- dat[, {
  .SD[bins, .(count= .N), .EACHI, on= c("seqnames", "start>=start", "start<=end")]
}, .(cdition, total, class, alt_class)]
setnames(counts, make.unique(names(counts)))
setnames(counts, "start.1", "end")
counts <- merge(counts[cdition!="PH18_2", !"end"],
                counts[cdition=="PH18_2", !"cdition"],
                by= c("class", "alt_class", "seqnames", "start"),
                suffixes= c("_mut", "_ctl"))
counts[, c("OR", "pval"):= {
  fisher.test(matrix(c(count_mut,
                       total_mut-count_mut,
                       count_ctl, 
                       total_ctl-count_ctl), 2, 2), 
              alternative='greater')[c("estimate", "p.value")]
}, .(count_mut, total_mut, count_ctl, total_ctl)]
counts[, fdr:= p.adjust(pval, "fdr"), .(cdition, alt_class)]
final <- vl_collapseBed(counts[fdr<0.1])

# Overlapping genes
genes <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)[type=="gene"]
ov <- vl_intersectBed(genes, final)[, .(gene_id, gene_symbol, seqnames, start, end)]
fwrite(ov,
       "Rdata/genes_overlapping_highly_mutated_regions.txt",
       sep= "\t",
       quote= F,
       na= NA)

# Plot
beds <- list.files("db/bed/mutations/", 
                   full.names = T)

pdf("pdf/gDNA_hostspots.pdf", 
    20, 
    4.5)
par(mar= c(5,10,5,2))
final[, {
  vl_screenshot(.SD,
                beds,
                genome= "dm6")
  title(main= paste0(.SD, collapse = "_"))
}, idx]
dev.off()

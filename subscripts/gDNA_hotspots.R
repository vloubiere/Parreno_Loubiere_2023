setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- fread("Rdata/gDNA_final_table.txt")
# dat <- dat[!(PH18)]
dat[, c("seqnames", "start"):= tstrsplit(id, ":|-", keep= c(1,2))]
dat[, seqnames:= paste0("chr", seqnames)]
dat[, start:= as.numeric(start)]
dat[, total:= .N, .(cdition, class, alt_class)]

bins <- vl_binBSgenome("dm6", 
                       bins_width = 25e3, 
                       restrict_seqnames = c("chr2L", "chr2R", "chr3L", "chr3R", "chr4","chrX","chrY"))

counts <- dat[, {
  .SD[bins, .(count= .N), .EACHI, on= c("seqnames", "start>=start", "start<=end")]
}, .(cdition, total, class, alt_class)]
setnames(counts, make.unique(names(counts)))
counts$start.1 <- NULL
counts <- merge(counts[cdition!="PH18_2"],
                counts[cdition=="PH18_2", !"cdition"],
                by= c("class", "alt_class", "seqnames", "start"),
                suffixes= c("_mut", "_ctl"))
counts[, pval:= {
  fisher.test(matrix(c(count_mut,
                       total_mut-count_mut,
                       count_ctl, 
                       total_ctl-count_ctl), 2, 2), alternative='greater')$p.value
}, .(count_mut, total_mut, count_ctl, total_ctl)]
counts[, fdr:= p.adjust(pval, "fdr"), .(cdition, alt_class)]
counts[fdr<0.1]

pdf("pdf/Figures/gDNA_hostspots.pdf", 
    15, 
    4.5)
par(cex= 0.6)
counts[, {
  ids <- split(id, cdition)
  vl_upset_plot(ids)
  title(main= paste(class, alt_class))
}, .(class, alt_class)]
dev.off()

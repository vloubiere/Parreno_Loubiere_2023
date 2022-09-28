setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)

# Import
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")
dat <- meta[!is.na(FC_peaks), fread(FC_peaks), .(ChIP, cdition, FC_peaks)]
dat[padj<0.05 & abs(log2FoldChange)>1, change:= paste0(ChIP, "_", cdition, "_",  ifelse(log2FoldChange>0, "up", "down")), .(cdition, ChIP)]
dat <- dat[grepl("H3K27me3.*down|H3K27Ac.*up|H3K4me1.*up", change)]
dat <- vl_importBed(dat[, c("seqnames", "start", "end"):= tstrsplit(ID, "_|:|-", keep= 2:4)])

coll <- vl_collapseBed(dat)
dat[coll, idx:= i.idx, on= c("seqnames", "start<=end", "end>=start")]
categ <- split(dat$idx, dat$change)

pdf("pdf/cutnrun_intersect_diff_regions.pdf", 10, 7)
par(mfrow= c(2,2),
    mar= c(7,18,2,2))
vl_upset_plot(categ[grep("PH29", names(categ), value = T)])
vl_upset_plot(categ[grep("H3K27me3", names(categ), value = T)])
vl_upset_plot(categ[grep("H3K27Ac", names(categ), value = T)])
vl_upset_plot(categ[grep("H3K4me1", names(categ), value = T)])
dev.off()

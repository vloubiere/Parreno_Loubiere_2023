setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
require(pheatmap)
require(rtracklayer)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Import TSSs
tss <- import("../../genomes/flybase/dmel-all-r6.35_simplified.gtf")
tss <- resize(tss[tss$type=="mRNA"], 1, "start")
tss <- tss[order(as.character(seqnames(tss)), start(tss))]
tmp1 <- tempfile(fileext = ".bed")
export(granges(tss), tmp1)
  
# Import ATAC signal
ATAC <- import.bw("db/bw/ATAC/FRT82_rep1_uniq.bw")
ATAC <- ATAC[as.character(seqnames(ATAC)) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")]
ATAC$name <- paste0("bin", seq(length(ATAC)))
ATAC <- ATAC[order(as.character(seqlevels(ATAC)), start(ATAC))]
tmp2 <- tempfile(fileext = ".bed")
export(ATAC, tmp2)

# Distance to tss
close <- fread(cmd= closestBed(tmp2, tmp1))

# signal
candidates <- close[V13>1000, .(seqnames= V1, start=V2, end=V3, name=V4, score= V5)]
candidates[, score:= score/(end-start)]
candidates <- GRanges(candidates[score>5])

# Merge peaks
peaks <- reduce(candidates, min.gapwidth= 500)
peaks <- resize(peaks, 1000, "center")
peaks <- peaks[start(peaks)>0]
peaks$name <- paste0("peak_", seq(length(peaks)))

# quantify K4me1 and K4me3 signals
ChIP <- data.table(file= list.files("db/bed/ChIP", "H3K4me1|H3K4me3|INPUTh_", full.names = T))
ChIP[, cdition:= tstrsplit(basename(file), "_", keep=1)]
ChIP[, cdition:= gsub("INPUTh", "INPUT", cdition)]
quantif <- ChIP[,my_countReads(peaks, file), (ChIP)]
res <- quantif[, .(norm_counts= (sum(counts)+1)/(sum(total_reads))*1e6), .(seqnames, start, end, name, cdition)]
res <- res[res[, .(cdition=="INPUT", name, norm_counts)], input_norm_counts := i.norm_counts, on= "name"]
res <- res[cdition!="INPUT"]
res[, log2_enr:= log2(norm_counts)-log2(input_norm_counts), name]
dmat <- dcast(res, seqnames+start+end+name~cdition, value.var = "log2_enr")
final <- dmat[H3K4me1>1 & H3K4me1-H3K4me3>1]
final[, name:= paste0("enhancer_", seq(nrow(final)))]
saveRDS(GRanges(final[, seqnames:name]), "Rdata/ED_enhancers.rds")





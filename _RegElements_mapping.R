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
tss$name <- tss$gene_id
tss$score <- NULL
tmp1 <- tempfile(fileext = ".bed")
export(tss, tmp1)

# Import ATAC signal
ATAC <- import.bw("db/bw/ATAC/FRT82_rep1_uniq.bw")
ATAC <- ATAC[as.character(seqnames(ATAC)) %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")]
ATAC$name <- paste0("bin", seq(length(ATAC)))
ATAC <- ATAC[order(as.character(seqnames(ATAC)), start(ATAC))]
tmp2 <- tempfile(fileext = ".bed")
export(ATAC, tmp2)

# signal
candidates <- as.data.table(copy(ATAC))
candidates[, score:= score/(end-start)]
candidates <- GRanges(candidates[score>5])

# Merge peaks
peaks <- reduce(candidates, min.gapwidth= 500)
peaks <- resize(peaks, 1000, "center")
peaks <- peaks[start(peaks)>0]
peaks$name <- paste0("peak_", seq(length(peaks)))

# quantify K4me1 and K4me3 enrichments
ChIP <- data.table(file= list.files("db/bed", "FRT82|H3K4me1|H3K4me3|INPUTh", recursive = T, full.names = T))
ChIP[, cdition:= tstrsplit(basename(file), "_", keep=1)]
quantif <- ChIP[,my_countReads(peaks, file), (ChIP)]
res <- quantif[, .(norm_counts= (sum(counts)+1)/(sum(total_reads))*1e6), .(seqnames, start, end, name, cdition)]
res <- res[res[, .(cdition=="INPUTh", name, norm_counts)], input_norm_counts := i.norm_counts, on= "name"][]
res[grepl("FRT82", cdition), median_counts:= median(norm_counts)]
res[grepl("^H3K4", cdition), log2_enr:= log2(norm_counts)-log2(input_norm_counts)]
res[grepl("FRT82", cdition), log2_enr:= log2(norm_counts)-log2(median_counts)]
res <- res[cdition!="INPUTh"]
dmat <- dcast(res, seqnames+start+end+name~cdition, value.var = "log2_enr")

# Sort regulatory elements
RE <- dmat[FRT82>1.5 | H3K4me1>1 | H3K4me3>1]
colnames(RE)[5] <- "ATAC"
RE[, K4me1_K4me3_ratio:= H3K4me1-H3K4me3]
RE[, name:= paste0("RE_", .I)]
saveRDS(resize(GRanges(RE), 1, "center"), "Rdata/RE_ED.rds")
export(resize(GRanges(RE), 1, "center"), "Rdata/RE_ED.bed")





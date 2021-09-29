genes <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
genes[, seqnames:= paste0("chr", seqnames)]
mRNA <- genes[type=="mRNA"]

# Promoter elements
prom <- copy(mRNA[strand %in% c("+", "-")]) # kick out few genes for which strand= "*" 
prom[strand=="+", c("start", "end"):= .(start-250, start+250)]
prom[strand=="-", c("start", "end"):= .(end-250, end+250)]
prom <- unique(prom[, .(seqnames, start, end, gene_id)])

# Non-promoter-overlapping REs
ATAC <- fread("../../public_data/dm6/peaks/merged_ATAC_seq_peaks.txt")
ATAC$ov_prom <- prom[ATAC, .N>0, .EACHI, on= c("seqnames==Chr", "end>=Start", "start<=End")]$V1
enh <- unique(ATAC[!(ov_prom), .(seqnames= Chr, start= Start, end= End)])
enh[, c("start", "end"):= .(round((start+end)/2)-250, round((start+end)/2)+250)]
# Assign closest promoter
enh$gene_id <- prom[enh, gene_id[which.min(abs((start+end)/2-(i.start+i.end)/2))], .EACHI, on= "seqnames"]$V1
# Add overlapping gene(s)
setkeyv(enh, c("seqnames", "start", "end"))
setkeyv(mRNA, c("seqnames", "start", "end"))
ov <- na.omit(unique(foverlaps(enh, mRNA)[, .(seqnames, start= i.start, end= i.end, gene_id)]))
enh <- unique(rbind(enh, ov))

# Merge promoters and REs + retieve symbols and strand
prom[, ID:= paste0("TSS", .SD[,.I]), gene_id]
enh[, ID:= paste0("RE", .SD[,.I]), gene_id]
final <- unique(rbind(prom, enh))
final <- final[seqnames %in% c("chrX", "chrY", "chr2L", "chr2R", "chr3L", "chr3R", "chr4")]
canonical <- genes[type=="gene", .(gene_id, 
                                   gene_symbol, 
                                   seqnames_gene= seqnames, 
                                   start_gene= start, 
                                   end_gene= end, 
                                   strand_gene= strand)]
final <- final[canonical, , on= "gene_id", nomatch= NULL]

# add PRC1 cluster from loubiere et al 2020
peaks <- fread("../../public_data/dm6/peaks/PcG_peaks_Loubiere_SA_2020.txt")
colnames(peaks)[1:3] <- c("seqnames", "start", "end")
setkeyv(peaks, c("seqnames", "start", "end"))
setkeyv(final, c("seqnames", "start", "end"))
ov <- foverlaps(peaks, final)
final[ov, PRC1_binding:= i.Cluster, on= c("gene_id", "ID")]
final <- final[, .(FBgn= gene_id,
                   symbol= gene_symbol,
                   gene_coor= paste0(seqnames_gene, ":", start_gene, "-", end_gene, ":", strand_gene),
                   RE_ID= ID,
                   PRC1_cluster= PRC1_binding,
                   RE_coor= paste0(seqnames, ":", start, "-", end))]

# SAVE
saveRDS(final, "Rdata/gene_REs.rds")

# Plot example
source("D:/_R_data/functions/my_screenshot_new.R")
pdf("pdf/Screenshot_PcG_biding_assignment_strat.pdf")
par(mar= c(2,7,2,2))
my_screenshot(bw_GR_list = list(PH= "../../public_data/dm6/bw/PH_ED_merge.bw", 
                                H3K4me1= "../../public_data/dm6/bw/H3K4me1_ED_merge.bw",
                                "ATAC-Seq"= "../../public_data/dm6/bw/ATAC_JJ_ED_merge.bw",
                                PRC1_peaks= GRanges(peaks),
                                all_Notch_promoters= GRanges(final[grepl("^TSS", RE_ID) & symbol=="N", RE_coor]),
                                all_Notch_enhancers= GRanges(final[grepl("^RE", RE_ID) & symbol=="N", RE_coor]),
                                bound_Notch_promoters= GRanges(final[grepl("^TSS", RE_ID) & !is.na(PRC1_cluster) & symbol=="N", RE_coor]),
                                bound_Notch_enhancers= GRanges(final[grepl("^RE", RE_ID) & !is.na(PRC1_cluster) & symbol=="N", RE_coor])),
              bed = GRanges("chrX", IRanges(3134870-20000, 3172221+20000)), 
              genome = "dm6")
dev.off()







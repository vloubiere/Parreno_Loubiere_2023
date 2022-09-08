setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)
require(BSgenome)

#############################
# Import promoters / enhancers and compute counts
#############################
# proms
gtf <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
gtf <- as.data.table(gtf)[type=="gene"]
gtf[, seqnames:= paste0("chr", seqnames)]
proms <- vl_resizeBed(gtf, "start", 200, 50)

# ATAC
enh <- vl_importBed("db/peaks/ATAC/ATAC_confident_peaks.narrowPeak", extraCols= "narrowPeak")
enh <- enh[signalValue>2 & qValue>20]
enh[, start:= start+peak]
enh[, end:= start]
TSS <- vl_resizeBed(gtf, "start", 0, 0)
cl <- vl_closestBed(enh, TSS)
enh[cl, gene_id:= i.gene_id.b, on= "name", mult= "first"]
enh <- vl_resizeBed(enh, "start", 125, 125)

# PH peask from SA2020
PH <- vl_importBed("external_data/PRC1_summits_SA2020_aax4001_table_s3.txt")
cl <- vl_closestBed(PH, TSS)
PH[cl, gene_id:= i.gene_id.b, on= "name", mult= "first"]
PH <- vl_resizeBed(PH, "start", 125, 125)

#############################
# Merge
#############################
dat <- rbindlist(list(TSS= proms[, .(gene_id,seqnames,start,end)],
                      ATAC= enh,
                      PH= PH),
                 idcol = "group",
                 fill= T)
dat <- dat[, .(seqnames, start, end, group, FBgn= gene_id)]
dat[gtf, symbol:= i.gene_symbol, on="FBgn==gene_id"]
dat <- dat[seqnames %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")]

#############################
# Motif counts
#############################
dat[, seq:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6, 
                                names= seqnames,
                                start= start,
                                end= end))]
counts <- vl_motif_counts(sequences = dat$seq)
counts <- as.data.table(counts)
counts <- cbind(dat[, .(group, coor= paste0(seqnames, ":", start, "-", end), symbol, FBgn)],
                counts)
fwrite(counts,
       "Rdata/final_RE_motifs_table.txt", 
       sep= "\t",
       quote= F,
       na= NA)

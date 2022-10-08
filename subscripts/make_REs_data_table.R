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
proms <- vl_resizeBed(gtf, "start", 200, 50, genome = "dm6")
proms[, dist:=0]

# ATAC
enh <- vl_importBed("db/peaks/ATAC/ATAC_confident_peaks.narrowPeak")
enh <- enh[signalValue>2 & qValue>20]
enh[, start:= start+peak]
enh[, end:= start]
TSS <- vl_resizeBed(gtf, "start", 0, 0)
cl <- vl_closestBed(enh, TSS)
cl <- cl[abs(dist)<25000]
enh[cl, c("gene_id", "dist"):= .(i.gene_id.b, i.dist), on= "name", mult= "first"]
enh <- enh[!is.na(gene_id)]
enh <- vl_resizeBed(enh, "center", 125, 125, genome = "dm6")

# PH peaks PH18
PH <- vl_importBed("db/peaks/cutnrun/PH_PH18_confident_peaks.narrowPeak")
cl <- vl_closestBed(PH, TSS)
PH[cl, c("gene_id", "dist"):= .(i.gene_id.b, i.dist), on= "name", mult= "first"]
PH <- PH[!is.na(gene_id)]
PH <- vl_resizeBed(PH, "center", 125, 125, genome = "dm6")

#############################
# Merge
#############################
dat <- rbindlist(list(TSS= proms[, .(gene_id, seqnames, start, end, dist)],
                      ATAC= enh,
                      PH= PH),
                 idcol = "group",
                 fill= T)
dat <- dat[, .(seqnames, start, end, group, FBgn= gene_id, dist)]
dat[gtf, symbol:= i.gene_symbol, on="FBgn==gene_id"]
dat <- dat[seqnames %in% c("chr2L","chr2R","chr3L","chr3R","chr4","chrX")]

#############################
# Motif counts
#############################
dat[, seq:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6, 
                                names= seqnames,
                                start= start,
                                end= end))]
counts <- vl_motif_counts(sequences = dat$seq, 
                          sel = vl_Dmel_motifs_DB_full[collection %in% c("bergman", "flyfactorsurvey", "jaspar"), motif])
counts <- as.data.table(counts)
counts <- cbind(dat[, .(group, coor= paste0(seqnames, ":", start, "-", end), symbol, FBgn, dist)],
                counts)
fwrite(counts,
       "Rdata/final_RE_motifs_table.txt", 
       sep= "\t",
       quote= F,
       na= NA)

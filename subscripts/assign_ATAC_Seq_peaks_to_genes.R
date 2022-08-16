setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)
require(BSgenome)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#------------------------------------------------#
# Import promoters / enhancers and compute counts
#------------------------------------------------#
# proms
proms <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
proms <- as.data.table(proms)
proms[, seqnames:= paste0("chr", seqnames)]
proms <- proms[seqnames %in% c("chrX","chrY","chr2L","chr2R","chr3L","chr3R","chr4")]
proms <- proms[type=="gene", .(gene_id,
                               gene_symbol,
                               seqnames, 
                               start= ifelse(strand=="+", start, end), 
                               end= ifelse(strand=="+", start, end), 
                               strand)]
# enh
enh <- vl_importBed("db/peaks/ATAC/ATAC_filtered_merged_peaks.narrowPeak")
enh <- enh[V7>2 & V9>20]
enh[, center:= rowMeans(.SD), .SDcols= c("start", "end")]
enh$closest_prom <- proms[enh, x.gene_id[which.min(abs(x.start-i.center))], .EACHI, on= "seqnames"]$V1
enh[proms, closest_symbol:= i.gene_symbol, on= "closest_prom==gene_id"]

# ph peask from SA2020
PH <- vl_importBed("external_data/PRC1_summits_SA2020_aax4001_table_s3.txt")
PH$closest_prom <- proms[PH, x.gene_id[which.min(abs(x.start-i.start))], .EACHI, on= "seqnames"]$V1
PH[proms, closest_symbol:= i.gene_symbol, on= "closest_prom==gene_id"]

# Merge
dat <- rbindlist(list(Promoter= proms[, .(FBgn= gene_id,
                                          symbol= gene_symbol,
                                          seqnames,
                                          start,
                                          end)],
                      ATAC= enh[, .(FBgn= closest_prom,
                                    symbol= closest_symbol,
                                    seqnames,
                                    start,
                                    end)],
                      PH= PH[, .(FBgn= closest_prom,
                                 symbol= closest_symbol,
                                 seqnames,
                                 start,
                                 end)]), idcol = "group")
dat <- vl_resizeBed(dat,
                    center = "center", 
                    upstream = 250, 
                    downstream = 250)
# Motif counts
dat[, seq:= as.character(getSeq(BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6, 
                                names= seqnames,
                                start= start,
                                end= end))]
counts <- vl_motif_counts(sequences = dat$seq)
colnames(counts) <- paste0(colnames(counts), "_mot")
counts <- cbind(dat,
                as.data.table(counts))
saveRDS(counts, 
        "Rdata/motif_counts_ATAC_sites.rds")

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Import clusters and define recovering genes
dat <- readRDS("Rdata/clustering_RNA.rds")$data
dat <- dat[cl %in% c(2,5), .(FBgn, RECOVERY= ifelse(cl==2, F, T))]
# Add gene coor
genes <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[type=="gene", .(FBgn= gene_id, seqnames, start, end, strand, symbol= gene_symbol)]
genes <- vl_resizeBed(genes, "start", 5000, genes[,end-start])
dat <- genes[dat, on= "FBgn"]
# Keep only K27me3+, PRC1 bound genes
dat <- dat[vl_covBed(dat, "external_data/PRC1_summits_SA2020_aax4001_table_s3.txt")>0]
dat <- dat[vl_covBed(dat, "db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.bed")>0]

setcolorder(dat,
            c("FBgn", "symbol"))
# SAVE
fwrite(dat, 
       "Rdata/cl2_cl5_RECOVERY_genes.txt",
       na = NA, 
       sep= "\t")

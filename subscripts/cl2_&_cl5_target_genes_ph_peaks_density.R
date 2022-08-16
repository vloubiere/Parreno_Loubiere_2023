setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

dat <- fread("Rdata/cl2_cl5_RECOVERY_genes.txt")[, .(FBgn, symbol, RECOVERY)]

# Add prom coor
genes <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[type=="gene", .(FBgn= gene_id, seqnames, start, end, strand, symbol= gene_symbol)]
proms <- vl_resizeBed(genes, "start", 0, 0)
dat <- proms[dat, on= "FBgn"]

# Distance to 5 closest PH peaks
PH <- vl_importBed("external_data/PRC1_summits_SA2020_aax4001_table_s3.txt")

dist <- vl_closestBed(dat, PH, n = 5)
dist <- dist[, .(mean_dist= mean(log10(dist+1)), min_dist= log10(min(dist+1))), .(seqnames, start, end)]

dat <- dist[dat, on= c("seqnames", "start", "end")]

vl_boxplot(mean_dist~RECOVERY, dat, compute_pval= list(c(1,2)))
vl_boxplot(min_dist~RECOVERY, dat, compute_pval= list(c(1,2)))

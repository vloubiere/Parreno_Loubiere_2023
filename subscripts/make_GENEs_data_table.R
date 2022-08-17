setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)

##########################################################
# Import RNA-Seq data
##########################################################
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA" & FC_file!="NA"]
meta[, cdition:= gsub("RNA_", "", cdition)]
meta[, cdition:= factor(cdition, c("PH18", "PHD11", "PHD9", "PH29"))]
dat <- meta[, fread(FC_file), .(cdition, FC_file)]
dat <- dcast(dat, 
             FBgn~cdition, 
             value.var = list("log2FoldChange", "padj"))

##########################################################
# Add RNA clusters
##########################################################
som <- readRDS("Rdata/clustering_RNA.rds")
dat[data.table(FBgn= rownames(som$data[[1]]), cl= som$unit.classif), cl:= i.cl, on= "FBgn"]

##########################################################
# Add gene symbols
##########################################################
gtf <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)[type=="gene"]
gtf <- gtf[dat$FBgn, on= "gene_id"]
dat[, symbol:= gtf$gene_symbol]
dat[, col:= vl_palette_few_categ(max(cl, na.rm = T))[cl]]

##########################################################
# PRC1 and K27me3 binding
##########################################################
# PRC1
PRC1 <- vl_importBed("external_data/PRC1_summits_SA2020_aax4001_table_s3.txt")
dat[, PRC1_bound:= vl_covBed(vl_resizeBed(gtf, "start", upstream = 2500, downstream = gtf$width+1), PRC1)>0]
# K27me3
K27me3 <- vl_importBed("db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.bed")
dat[, K27me3_bound:= vl_covBed(vl_resizeBed(gtf, "start", upstream = 2500, downstream = gtf$width+1), K27me3)>0]

##########################################################
# Add RECOVERY
##########################################################
dat[cl==2 & PRC1_bound & K27me3_bound, recovery:= F]
dat[cl==5 & PRC1_bound & K27me3_bound, recovery:= T]

##########################################################
# Quantif CUTNRUN
##########################################################
# Add HTMs CUTNRUN and CHIP-Seq
files <- list.files("db/bw/cutnrun/", "_merge.bw", full.names = T)
files <- c(files,
           list.files("db/bw/SA_2020/", "PC_ED|PH_ED|SUZ12_ED", full.names = T))
.n <- gsub("_merge.bw", "", basename(files))
TSSs <- vl_resizeBed(gtf, "start", upstream = 500, downstream = 500)
dat[, paste0(.n, "_TSS"):= lapply(files, function(x) vl_bw_coverage(TSSs, x))]
BODYs <- vl_resizeBed(gtf, "start", upstream = 2500, downstream = gtf$width+1)
dat[, paste0(.n, "_BODY"):= lapply(files, function(x) vl_bw_coverage(BODYs, x))]

##########################################################
# Add chromatin types SA2020
##########################################################
chromhmm <- vl_importBed("external_data/chromatin_types_SA2020_table_s1.txt")
gtf[chromhmm, chromhmm:= i.name, on= c("seqnames", "start<=end", "start>=start"), mult= "first"]
dat$chromhmm <- gtf$chromhmm

setcolorder(dat, 
            c("FBgn", "symbol", "PRC1_bound", "K27me3_bound", "cl", "recovery", "chromhmm", "col"))
fwrite(dat,
       "Rdata/final_gene_features_table.txt", 
       sep= "\t",
       quote= F,
       na= NA)


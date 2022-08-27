setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)

##########################################################
# Import RNA-Seq data
##########################################################
meta <- fread("Rdata/processed_metadata_RNA_align.txt")
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
som <- readRDS("Rdata/clustering_RNA_align.rds")
dat[data.table(FBgn= rownames(som$data[[1]]), cl= som$unit.classif), cl:= i.cl, on= "FBgn"]

##########################################################
# Add FPKMs
##########################################################
dds <- readRDS("db/dds_align/epiCancer_ED_GFP-_system_RNA_dds.rds")
fpkms <- as.data.table(DESeq2::fpkm(dds), keep.rownames = "FBgn")
fpkms <- melt(fpkms, id.vars = "FBgn")
fpkms[, variable:= tstrsplit(variable, "_", keep=1)]
fpkms <- fpkms[, .(FPKM= mean(value)), .(variable, FBgn)]
fpkms <- dcast(fpkms, FBgn~variable, value.var = "FPKM")
setnames(fpkms, 
         c("PH18", "PH29", "PHJ11", "PHJ9", "W18", "W29", "WKD"), 
         c("FPKM_PH18", "FPKM_PH29", "FPKM_PHD11", "FPKM_PHD9", "FPKM_W18", "FPKM_W29", "FPKM_WKD"))
dat <- fpkms[dat, on= "FBgn"]

##########################################################
# Add gene symbols
##########################################################
gtf <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
gtf <- as.data.table(gtf)[type=="gene", .(FBgn= gene_id, symbol= gene_symbol, seqnames, start, end, strand)]
dat <- gtf[dat, on= "FBgn"]
dat[, col:= vl_palette_few_categ(max(cl, na.rm = T))[cl]]

##########################################################
# PRC1 and K27me3 binding
##########################################################
PRC1 <- unlist(get(load("external_data/SA2020_cl.list"))$genes)
dat[, PRC1_bound:= symbol %in% PRC1]
# K27me3
K27me3 <- vl_importBed("db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak", specialFormat = "broadPeak")
dat[, K27me3_bound:= vl_covBed(vl_resizeBed(data.table(seqnames, start, end, strand), 
                                            "start", 
                                            upstream = 2500,
                                            downstream = end-start+1), K27me3)>0]

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
TSSs <- vl_resizeBed(dat[, .(seqnames, start, end)], 
                     "start", 
                     upstream = 500, 
                     downstream = 500)
dat[, paste0(.n, "_TSS"):= lapply(files, function(x) vl_bw_coverage(TSSs, x))]
BODYs <- vl_resizeBed(dat[, .(seqnames, start, end)], 
                      "start", 
                      upstream = 2500, 
                      downstream = dat[, end-start+1])
dat[, paste0(.n, "_BODY"):= lapply(files, function(x) vl_bw_coverage(BODYs, x))]

##########################################################
# Add chromatin types SA2020
##########################################################
chromhmm <- vl_importBed("external_data/chromatin_types_SA2020_table_s1.txt")
TSSs <- vl_resizeBed(dat[, .(seqnames, start, end)], 
                     "start", 
                     upstream = 0, 
                     downstream = 0)
TSSs[chromhmm, chromhmm:= i.name, on= c("seqnames", "start<=end", "start>=start"), mult= "first"]
dat$chromhmm <- TSSs$chromhmm

setcolorder(dat,
            c("FBgn", "symbol", "PRC1_bound", "K27me3_bound", "cl", "recovery", "chromhmm", "col", "seqnames", "start", "end", "strand"))
fwrite(dat,
       "Rdata/final_gene_features_table_align.txt", 
       sep= "\t",
       quote= F,
       na= NA)


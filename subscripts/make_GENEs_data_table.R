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
meta <- meta[DESeq2_object=="epiCancer_noGFP" & FC_file!="NA"]
meta[, cdition:= factor(cdition, c("PH18", "PHD11", "PHD9", "PH29"))]
dat <- meta[, fread(FC_file), .(cdition, FC_file)]
dat <- dcast(dat, 
             FBgn~cdition, 
             value.var = list("log2FoldChange", "padj", "diff"))

##########################################################
# Add RNA clusters
##########################################################
som <- readRDS("Rdata/clustering_RNA.rds")
dat[data.table(FBgn= rownames(som$data[[1]]), cl= som$unit.classif), cl:= i.cl, on= "FBgn"]
dat[, col:= vl_palette_few_categ(max(cl, na.rm = T))[cl]]

##########################################################
# Add FPKMs
##########################################################
dds <- readRDS("db/dds/RNA/epiCancer_noGFP_dds.rds")
fpkms <- as.data.table(DESeq2::fpkm(dds), keep.rownames = "FBgn")
fpkms <- melt(fpkms, id.vars = "FBgn")
fpkms[, variable:= tstrsplit(variable, "_", keep=1)]
fpkms <- fpkms[, .(FPKM= mean(value)), .(variable, FBgn)]
fpkms <- dcast(fpkms, FBgn~variable, value.var = "FPKM")
setnames(fpkms, names(fpkms)[-1], paste0("FPKM_", names(fpkms)[-1]))
dat <- fpkms[dat, on= "FBgn"]

##########################################################
# Add gene coordinates and symbols
##########################################################
gtf <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
dat[gtf[type=="gene"], c("symbol", "seqnames", "start", "end", "strand"):= 
      .(i.gene_symbol, seqnames, start, end, strand), on= "FBgn==gene_id"]

##########################################################
# PRC1 and K27me3 binding
##########################################################
PRC1 <- unlist(get(load("external_data/SA2020_cl.list"))$genes)
dat[, PRC1_bound:= symbol %in% PRC1]
# K27me3
K27me3 <- vl_importBed("db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak", 
                       extraCols = "broadPeak")
K27me3 <- K27me3[signalValue>2]
dat[, K27me3_bound:= vl_covBed(dat, K27me3)>0]

##########################################################
# Add RECOVERY
##########################################################
dat[, unaffected_PH18:= (is.na(padj_PH18) | padj_PH18>0.05)]
dat[, up_PH29:= padj_PH29<0.05 & log2FoldChange_PH29>(log2(1.5))]
dat[, up_PHD9:= padj_PHD9<0.05 & log2FoldChange_PHD9>(log2(1.5))]
dat[, up_PHD11:= padj_PHD11<0.05 & log2FoldChange_PHD11>(log2(1.5))]
dat[K27me3_bound
    & unaffected_PH18
    & up_PH29, recovery:= ifelse(up_PHD9 | up_PHD11, FALSE, TRUE)]
dat$unaffected_PH18 <- dat$up_PH29 <- dat$up_PHD9 <- dat$up_PHD11 <- NULL

##########################################################
# Quantif CUTNRUN
##########################################################
# Body marks
gene_body <- vl_resizeBed(dat, 
                          center = "start", 
                          upstream = 1000, 
                          downstream = dat[, end-start+1], 
                          genome = "dm6")
files <- list.files("db/bw/cutnrun/", 
                    "^H2AK118Ub.*merge|^H3K27me3.*merge|^H3K36me3.*merge|^H3K4me1.*merge", 
                    full.names = T)
.n <- paste0(gsub("_merge.bw$", "", basename(files)), "_body")
dat[, (.n):= lapply(files, function(x) vl_bw_coverage(gene_body, x))]

# Promoter marks
prom <- vl_resizeBed(dat, 
                     center = "start", 
                     upstream = 750, 
                     downstream = 250, 
                     genome = "dm6")
files <- list.files("db/bw/SA_2020/", "PC_ED|PH_ED|SUZ12_ED", full.names = T)
files <- c(files,
           list.files("db/bw/cutnrun/", "^PH.*merge|^H3K27Ac.*merge", full.names = T))
.n <- paste0(gsub("_merge.bw$", "", basename(files)), "_prom")
dat[, (.n):= lapply(files, function(x) vl_bw_coverage(prom, x))]

##########################################################
# Add chromatin types SA2020
##########################################################
chromhmm <- vl_importBed("external_data/chromatin_types_SA2020_table_s1.txt")
dat[chromhmm, chromhmm:= i.name, .EACHI, on= c("seqnames", "start<=end", "start>=start")]


##########################################################
# Save
##########################################################
setcolorder(dat,
            c("FBgn", "symbol", "PRC1_bound", "K27me3_bound", "cl", 
              "recovery", "chromhmm", "col", "seqnames", "start", "end", "strand", 
              "diff_PH18", "diff_PH29", "diff_PHD9", "diff_PHD11"))
fwrite(dat,
       "Rdata/final_gene_features_table.txt", 
       sep= "\t",
       quote= F,
       na= NA)

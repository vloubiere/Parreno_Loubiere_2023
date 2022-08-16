setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

############################################################
# Make tables AMM
############################################################
# Transcriptomes ------------------------------------------#
meta <- fread("Rdata/processed_metadata_RNA.txt")
files <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA", c(FC_file, FC_file_PH29)]
files <- unique(na.omit(files))
tab <- lapply(files, fread)
names(tab) <- gsub("epiCancer_ED_|RNA_|.txt", "", basename(files))
tab <- rbindlist(tab, idcol = "cdition")
tab <- dcast(tab,
             FBgn~cdition, 
             value.var = list("log2FoldChange", "padj"))
# FPKMs ---------------------------------------------------#
dds <- readRDS("db/dds/RNA/epiCancer_ED_GFP-_system_RNA_dds.rds")
FPKM <- as.data.table(DESeq2::fpkm(dds), keep.rownames = "FBgn")
FPKM <- melt(FPKM, id.vars= "FBgn")
FPKM[, variable:= tstrsplit(variable, "_", keep = 1)]
FPKM[, variable:= paste0(variable, "_FPKM")]
FPKM <- FPKM[, .(value= mean(value)), .(FBgn, variable)]
FPKM <- dcast(FPKM, FBgn~variable, value.var = "value")
tab<- tab[FPKM, on="FBgn"]
# Add symbols ---------------------------------------------#
gtf <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
gtf <- gtf[type=="gene", .(FBgn= gene_id, symbol= gene_symbol, seqnames, start= ifelse(strand=="-", end, start), strand)]
gtf[, end:= start]
tab <- gtf[tab, on= "FBgn"]
# Add clusters --------------------------------------------#
cl <- readRDS("Rdata/clustering_RNA.rds")
tab[cl$data, RNA_cluster:= i.cl, on="FBgn"]
recov <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")[, .(FBgn, PHD9_RECOVERY, PHD11_RECOVERY)]
tab <- recov[tab, on="FBgn"]
# Add PRC1 binding ----------------------------------------#
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
PRC1 <- loadRData("external_data/SA2020_cl.list")
PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
tab[PRC1, PRC1_cluster:= i.PRC1_cluster, on= "symbol"]
# HTMs binding --------------------------------------------#
HTM <- list.files("db/peaks/cutnrun/", "confident_peaks.bed", full.names = T)
names(HTM) <- gsub("_confident_peaks.bed", "", basename(HTM))
HTM <- rbindlist(lapply(HTM, fread, sel= 1:3, col.names= c("seqnames", "start", "end")), idcol = "cdition")
HTM[, {
  .c <- copy(.SD)
  bound <- tab[.c, FBgn, on= c("seqnames", "start<=end", "end>=start")]
  tab[, (cdition):= FBgn %in% bound]
  print(".")
}, cdition]

############################################################
# ORDER and SAVE
############################################################
setnames(tab, gsub("log2FoldChange", "log2FC", names(tab)))
setcolorder(tab, 
            c("FBgn", "symbol", "RNA_cluster", "PRC1_cluster", "PHD9_RECOVERY", "PHD11_RECOVERY"))
fwrite(tab, 
       "db/FC_tables/RNA_table_AMM.txt",
       col.names = T, 
       row.names = F, 
       sep= "\t",
       quote= F, 
       na = NA)

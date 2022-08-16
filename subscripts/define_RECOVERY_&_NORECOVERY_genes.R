setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import genes which are up-regulated at 29
dat <- data.table(file= c("db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PH29_vs_RNA_W29.txt",
                          "db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD9_vs_RNA_PH29.txt",
                          "db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD11_vs_RNA_PH29.txt",
                          "db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD9_vs_RNA_WKD.txt",
                          "db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD11_vs_RNA_WKD.txt"),
                  cdition= c("PH29_vs_W29",
                             "PHD9_vs_PH29",
                             "PHD11_vs_PH29",
                             "PHD9_vs_WKD",
                             "PHD11_vs_WKD"))
dat <- dat[, fread(file, sel= c(1,3,6)), cdition]
dat <- dcast(dat,
             FBgn~cdition,
             value.var = list("log2FoldChange", "padj"))
dat <- dat[padj_PH29_vs_W29<0.05 & log2FoldChange_PH29_vs_W29>log2(1.5)]
# Add gene coor
proms <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(proms) <- "UCSC"
proms <- as.data.table(proms)
proms <- proms[type=="gene", .(FBgn= gene_id, seqnames, start, end, strand, symbol= gene_symbol)]
proms <- vl_resizeBed(proms, "start", 2500, 2500)
dat <- proms[dat, on= "FBgn"]
# Keep only PRC1 bound genes
PRC1 <- loadRData("external_data/SA2020_cl.list")
PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
dat <- dat[symbol %in% PRC1$symbol]
# Keep only K27me3 genes
dat <- dat[vl_covBed(dat, "db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.bed")>0]
# Define recovering and non recovering sets of genes
dat[, PHD9_RECOVERY:= fcase(log2FoldChange_PHD9_vs_PH29<(-1) & padj_PHD9_vs_PH29<0.05, T,
                            log2FoldChange_PHD9_vs_WKD>0, F,
                            default = NA)]
dat[, PHD11_RECOVERY:= fcase(log2FoldChange_PHD11_vs_PH29<(-1) & padj_PHD11_vs_PH29<0.05, T,
                             log2FoldChange_PHD11_vs_WKD>0, F,
                             default = NA)]
setcolorder(dat,
            c("FBgn", "symbol"))
# SAVE
fwrite(dat, 
       "Rdata/RECOVERY_NORECOVERY_genes.txt",
       na = NA, 
       sep= "\t")

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

###################################################
# Format dat
###################################################
# Import
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA"]
meta <- na.omit(meta[, .(cdition= gsub("^RNA_", "", cdition), FC_file)])
meta[, cdition:= factor(cdition, 
                        levels= c("PH18", "PH29", "PHD9", "PHD11"))]
dat <- meta[, fread(FC_file), (meta)]
dat <- dat[diff!="unaffected"]
symbols <- as.data.table(rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf"))
symbols <- unique(symbols[, .(FBgn= gene_id, symbol= gene_symbol)])
dat <- symbols[dat, on= "FBgn", nomatch= NULL]
# Select diff genes no diff in PH18
dat <- dat[diff!="unaffected"]
diff_18 <- dat["PH18", FBgn, on= "cdition"]
dat <- dat[!(FBgn %in% diff_18) & cdition!="PH18"]
# Spli PRC1 +/-
PRC1 <- loadRData("external_data/SA2020_cl.list")
PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
dat[, class:= paste(diff, ifelse(symbol %in% PRC1$symbol, "PRC1+", "PRC1-"))]

###################################################
# plot
###################################################
pdf("pdf/Figures/GO_up_down_genes_GFP-.pdf", 
    width = 10, 
    height = 10)
par(las= 2)
dat[, {
  .c <- vl_GO_enrich(geneIDs = split(FBgn, class), 
                     species = "Dm", 
                     plot = F)
  plot(.c,
       padj_cutoff = 0.05,
       top_enrich = 10,
       cex.balloons= 0.5, 
       main= cdition)
  print("done")
}, cdition]
dev.off()

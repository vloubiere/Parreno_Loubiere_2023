setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(scipen= 99)
require(pheatmap)
require(rtracklayer)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

#----------------------------------------------#
# 1- ChIP enrichment
#----------------------------------------------#
gtf <- as.data.table(import("/groups/stark/vloubiere/genomes/flybase/dmel-all-r6.35_simplified.gtf"))
# Gene body
g <- gtf[type=="mRNA", 
         .(seqnames, start, end, strand, name= paste(c(gene_id, gene_symbol, transcript_id, transcript_symbol), collapse= "__")), 
         .(gene_id, gene_symbol, transcript_id, transcript_symbol)]
g <- GRanges(g)
# tss
tss <- resize(g, 1, "start")
tss <- resize(tss, 500, "center")
tss <- tss[order(as.character(seqnames(tss)), start(tss))]
tmp1 <- tempfile(pattern = "tss", fileext = ".bed")
export(tss, tmp1)
# Closest enhancer
enh <- readRDS("Rdata/ED_enhancers.rds")
enh <- resize(enh, 500, "center")
tmp2 <- tempfile(pattern = "enh", fileext = ".bed")
export(enh, tmp2)
close <- fread(cmd= closestBed(tmp1, tmp2))
enh <- GRanges(close[V7!=".", .(seqnames= V7[1], start= V8[1], end= V9[1]), .(name= V4)])
# Compile
g$name <- paste0(g$name, "__gene")
tss$name <- paste0(tss$name, "__tss")
enh$name <- paste0(enh$name, "__cenh")
genes <- c(g, tss, enh)
genes <- genes[order(seqnames(genes), start(genes))]
tmp <- tempfile(fileext = ".bed")
export(genes, tmp)

ChIP <- data.table(file= list.files("db/bed/ChIP", ".bed", full.names = T))
ChIP[, cdition:= gsub("_uniq.bed|_rep1|_rep2", "", basename(file))]
ChIP[cdition %in% c("PC", "PH", "H3K27me3"), input_group:= "INPUTa"]
ChIP[cdition %in% c("PSC", "SUZ12"), input_group:= "INPUTv"]
ChIP[cdition %in% c("EYGSE112868", "POLIIGSE112868"), input_group:= "INPUTGSE112868"]
ChIP[grepl("^H3|^H4|^H2", cdition) & is.na(input_group), input_group:= "INPUTh"]
if(!file.exists("Rdata/ChIP_quantif.rds"))
{
  quantif <- ChIP[, my_countReads(genes, file), c(colnames(ChIP))]
  saveRDS(quantif, "Rdata/ChIP_quantif_all.rds")
}
if(!file.exists("Rdata/ChIP_normalized_enrichment_tss_body.rds"))
{
  quantif <- readRDS("Rdata/ChIP_quantif_all.rds")
  norm <- quantif[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), .(name, cdition, input_group)]
  norm[norm, input_norm_counts:= i.norm_counts, on= c("name", "input_group==cdition")]
  norm <- norm[!grepl("INPUT", cdition)]
  norm[, enr:= log2(norm_counts)-log2(input_norm_counts)]
  norm[, c("FBgn", "symbol", "transcript_id", "transcript_symbol", "gene_part"):= tstrsplit(name, "__")]
  dat <- norm[, .(tss_enr= max(.SD[gene_part=="tss", enr]), 
                  gene_enr= mean(.SD[gene_part=="gene", enr]),
                  enh_enr= mean(.SD[gene_part=="cenh", enr])), .(cdition, FBgn)]
  dmat <- dcast(dat, FBgn~cdition, value.var = c("tss_enr", "gene_enr", "enh_enr"))
  saveRDS(dmat, "Rdata/ChIP_normalized_enrichment_tss_body_cenh.rds")
}


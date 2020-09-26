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
g <- gtf[type=="gene", .(seqnames, start, end, strand, name= paste0(gene_id, "__", gene_symbol))]
g <- g[order(seqnames, start)]
g <- GRanges(g)
# tss
tss <- gtf[type=="mRNA", 
           .(seqnames, start, end, strand, name= paste0(gene_id, "__", gene_symbol))]
tss <- tss[order(seqnames, start)]
tss <- resize(resize(GRanges(tss), 1, "start"), 1000, "center")
# Close enhancers
ext_g <- gtf[type=="gene", .(seqnames, start, end, strand, name= paste0(gene_id, "__", gene_symbol))]
ext_g[strand=="+", start:= ifelse(start-2000>1, start-2000, 1)]
ext_g[strand=="-", end:= end+2000]
ext_g <- GRanges(ext_g)
tmp1 <- tempfile(pattern = "ext_g", fileext = ".bed")
export(ext_g, tmp1)
enh <- readRDS("Rdata/ED_enhancers.rds")
tmp2 <- tempfile(pattern = "enh", fileext = ".bed")
export(enh, tmp2)
c_enh <- fread(cmd= paste("/software/2020/software/bedtools/2.27.1-foss-2018b/bin/intersectBed -wb -a",  tmp1, "-b", tmp2))
c_enh <- GRanges(c_enh[, .(seqnames= V7, start= V8, end= V9, name= V4)])
c_enh <- resize(c_enh, 1000, "center")
# Compile
g$name <- paste0(g$name, "__gene")
tss$name <- paste0(tss$name, "__tss")
c_enh$name <- paste0(c_enh$name, "__cenh")
genes <- c(g, tss, c_enh)
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
if(!file.exists("Rdata/ChIP_normalized_enrichment_tss_body_cenh.rds"))
{
  quantif <- readRDS("Rdata/ChIP_quantif_all.rds")
  norm <- quantif[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), .(name, cdition, input_group)]
  norm[norm, input_norm_counts:= i.norm_counts, on= c("name", "input_group==cdition")]
  norm <- norm[!grepl("INPUT", cdition)]
  norm[, enr:= log2(norm_counts)-log2(input_norm_counts)]
  norm[, c("FBgn", "symbol", "gene_part"):= tstrsplit(name, "__")]
  dat <- norm[, .(enr= log2(sum(2^enr))), .(cdition, FBgn, symbol, gene_part)]
  dmat <- dcast(dat, FBgn+symbol~cdition+gene_part, value.var = "enr")
  saveRDS(dmat, "Rdata/ChIP_normalized_enrichment_tss_body_cenh.rds")
}


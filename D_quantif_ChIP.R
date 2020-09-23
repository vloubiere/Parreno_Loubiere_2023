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
g <- gtf[type=="mRNA", 
         .(seqnames, start, end, strand, name= paste(c(gene_id, gene_symbol, transcript_id, transcript_symbol), collapse= "__")), 
         .(gene_id, gene_symbol, transcript_id, transcript_symbol)]
g <- GRanges(g)
tss <- resize(g, 1, "start")
tss <- resize(tss, 500, "center")
g$name <- paste0(g$name, "__gene")
tss$name <- paste0(tss$name, "__tss")
genes <- c(g, tss)
genes <- genes[order(seqnames(genes), start(genes))]
tmp <- tempfile(fileext = ".bed")
export(genes, tmp)

ChIP <- data.table(file= list.files("db/bed/ChIP/", ".bed", full.names = T))
ChIP[, cdition:= gsub("_uniq.bed|_rep1|_rep2", "", basename(file))]
ChIP[cdition %in% c("PC", "PH", "H3K27me3"), input_group:= "INPUTa"]
ChIP[cdition %in% c("PSC", "SUZ12"), input_group:= "INPUTv"]
ChIP[grepl("^H3|^H4|^H2", cdition) & is.na(input_group), input_group:= "INPUTh"]
if(!file.exists("Rdata/ChIP_quantif.rds"))
{
  quantif <- ChIP[, my_countReads(genes, file), c(colnames(ChIP))]
  saveRDS(quantif, "Rdata/ChIP_quantif.rds")
}
if(!file.exists("Rdata/ChIP_normalized_enrichment_tss_body.rds"))
{
  quantif <- readRDS("Rdata/ChIP_quantif.rds")
  norm <- quantif[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), .(name, cdition, input_group)]
  norm[norm, input_norm_counts:= i.norm_counts, on= c("name", "input_group==cdition")]
  norm <- norm[!grepl("INPUT", cdition)]
  norm[, enr:= log2(norm_counts)-log2(input_norm_counts)]
  norm[, c("FBgn", "symbol", "transcript_id", "transcript_symbol", "gene_part"):= tstrsplit(name, "__")]
  dat <- norm[, .(tss_enr= max(.SD[gene_part=="tss", enr]), gene_enr= mean(.SD[gene_part=="gene", enr])), .(cdition, FBgn)]
  dmat <- dcast(dat, FBgn~cdition, value.var = c("tss_enr", "gene_enr"))
  saveRDS(dmat, "Rdata/ChIP_normalized_enrichment_tss_body.rds")
}


setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(scipen= 99)
require(pheatmap)
require(rtracklayer)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

#----------------------------------------------#
# 1- Import peaks, import TSSs perform RE/gene assignment
#----------------------------------------------#
gtf <- as.data.table(import("/groups/stark/vloubiere/genomes/flybase/dmel-all-r6.35_simplified.gtf"))
# transcripts tss
tss <- gtf[type=="mRNA", .(seqnames, start, end, strand, gene_id, gene_symbol, transcript_id)]
tss <- resize(resize(GRanges(tss), 1, "start"), 1000, "center")
# canonical genes
cgenes <- gtf[type=="gene", .(seqnames, start, end, strand, name= gene_id)]
cgenes <- cgenes[order(seqnames, start)]
# Closest REs
RE <- readRDS("Rdata/RE_ED.rds")
RE <- RE[order(as.character(seqnames(RE)), start(RE))]
tmp1 <- tempfile(fileext = ".bed")
tmp2 <- tempfile(fileext = ".bed")
export(cgenes, tmp1)
export(RE, tmp2)
# Assign canonical genes with REs
ass <- fread(cmd= paste("/software/2020/software/bedtools/2.27.1-foss-2018b/bin/closestBed -D ref -a", tmp1, "-b", tmp2, "-k", 5))
ass <- ass[, .(seqnames= V7, start= V8, end= V9, FBgn= V4, dist= V13)]
ass <- ass[dist > -5000 & dist < 5000, !"dist"]
ass[, ID:= paste0("RE", .I), FBgn]
ass <- GRanges(ass[order(seqnames, start)])
# Compile
tss$name <- paste0(tss$gene_id, "_", tss$gene_symbol, "_", "TSS", seq(length(tss)))
ass$name <- paste0(ass$FBgn, "_", ass$ID)
genes <- c(tss, ass)

#----------------------------------------------#
# 2- Compute ChIP-Seq signal
#----------------------------------------------#
ChIP <- data.table(file= list.files("db/bed/ChIP/", ".bed", full.names = T))
ChIP[, cdition:= gsub("_uniq.bed|_rep1|_rep2", "", basename(file))]
ChIP[cdition %in% c("PC", "PH", "H3K27me3"), input_group:= "INPUTa"]
ChIP[cdition %in% c("PSC", "SUZ12"), input_group:= "INPUTv"]
ChIP[grepl("^POLII|^EY", cdition), input_group:= "INPUTGSE112868"]
ChIP[grepl("^H3|^H4|^H2", cdition) & is.na(input_group), input_group:= "INPUTh"]
if(!file.exists("Rdata/ChIP_quantif_TSSs_REs.rds"))
{
  quantif <- ChIP[, my_countReads(genes, file), c(colnames(ChIP))]
  saveRDS(quantif, "Rdata/ChIP_quantif_TSSs_REs.rds")
}
if(!file.exists("Rdata/ChIP_normalized_enrichment_TSS_REs.rds"))
{
  quantif <- readRDS("Rdata/ChIP_quantif_TSSs_REs.rds")
  norm <- quantif[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), .(name, cdition, input_group)]
  norm[norm, input_norm_counts:= i.norm_counts, on= c("name", "input_group==cdition")]
  norm <- norm[!grepl("INPUT", cdition)]
  norm[, enr:= log2(norm_counts)-log2(input_norm_counts)]
  norm[grepl("_TSS", name), c("FBgn", "symbol", "RE"):= tstrsplit(name, "_")]
  norm[!grepl("_TSS", name), c("FBgn", "RE"):= tstrsplit(name, "_")]
  tss <- norm[grepl("^TSS", RE), .(TSS_enr= max(enr, na.rm= T)), .(FBgn, symbol, cdition)]
  RE <- norm[, .(RE_enr= mean(enr, na.rm= T)), .(FBgn, cdition)]
  res <- merge(tss, RE, by= c("FBgn", "cdition"))
  dmat <- dcast(res, FBgn+symbol~cdition, value.var = c("TSS_enr", "RE_enr"))
  saveRDS(dmat, "Rdata/ChIP_normalized_enrichment_TSS_REs.rds")
}


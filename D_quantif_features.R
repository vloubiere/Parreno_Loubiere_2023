setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(datatable.print.topn= 1)
options(scipen= 99)
require(GenomicRanges)
require(rtracklayer)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
require(BSgenome.Dmelanogaster.UCSC.dm6)
require(org.Dm.eg.db)
require(pheatmap)
require(seqinr)
require(motifmatchr)
require(PWMEnrich)
require(TFBSTools)
require(seqLogo)
require(kohonen)

#----------------------------------------------#
# 1- Import peaks, import TSSs perform RE/gene assignment
#----------------------------------------------#
if(!file.exists("Rdata/genes_TSSs_REs.rds"))
{
  gtf <- as.data.table(import("/groups/stark/vloubiere/genomes/flybase/dmel-all-r6.35_simplified.gtf"))
  # transcripts tss
  tss <- gtf[type=="mRNA", .(seqnames, start, end, strand, gene_id, gene_symbol, transcript_id)]
  tss <- resize(resize(GRanges(tss), 1, "start"), 1000, "center")
  # canonical genes
  cgenes <- gtf[type=="gene" & gene_id %in% tss$gene_id, .(seqnames, start, end, strand, name= gene_id)]
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
  ass$name <- paste0(ass$FBgn, "_", tss$gene_symbol[match(ass$FBgn, tss$gene_id)], "_", ass$ID)
  genes <- c(tss, ass)
  genes <- resize(genes, 1000, "center")
  saveRDS(genes, "Rdata/genes_TSSs_REs.rds")
}

#----------------------------------------------#
# 2- Compute ChIP-Seq signal
#----------------------------------------------#
if(!file.exists("Rdata/ChIP_quantif_TSSs_REs.rds"))
{
  ChIP <- data.table(file= list.files("db/bed/ChIP/", ".bed", full.names = T))
  ChIP[, cdition:= gsub("_uniq.bed|_rep1|_rep2", "", basename(file))]
  ChIP[cdition %in% c("PC", "PH", "H3K27me3"), input_group:= "INPUTa"]
  ChIP[cdition %in% c("PSC", "SUZ12"), input_group:= "INPUTv"]
  ChIP[grepl("^POLII|^EY", cdition), input_group:= "INPUTGSE112868"]
  ChIP[grepl("^H3|^H4|^H2", cdition) & is.na(input_group), input_group:= "INPUTh"]
  genes <- readRDS("Rdata/genes_TSSs_REs.rds")
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
  tss <- norm[grepl("^TSS", RE), .(TSS_enr= max(enr, na.rm= T)), .(FBgn, symbol, cdition)]
  RE <- norm[, .(RE_enr= mean(enr, na.rm= T)), .(FBgn, cdition)]
  res <- merge(tss, RE, by= c("FBgn", "cdition"))
  dmat <- dcast(res, FBgn+symbol~cdition, value.var = c("TSS_enr", "RE_enr"))
  saveRDS(dmat, "Rdata/ChIP_normalized_enrichment_TSS_REs.rds")
}

#----------------------------------------------#
# 3- Collapse informative motifs
#----------------------------------------------#
if(!file.exists("Rdata/som_informative_motifs.rds"))
{
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  sel <- TF_clusters_PWMs$metadata[!is.na(TF_clusters_PWMs$metadata$Dmel), "motif_name"]
  sel <- match(sel, name(TF_clusters_PWMs$All_pwms_log_odds))
  # counts
  genes <- readRDS("Rdata/genes_TSSs_REs.rds")
  genes <- resize(genes, 250, "center")
  hit <- matchMotifs(TF_clusters_PWMs$All_pwms_log_odds[sel], genes, genome= "dm6",
                     p.cutoff= 5e-4, bg="even", out= "scores")
  counts <- as.matrix(motifCounts(hit))
  counts <- log2(counts+1)
  rownames(counts) <- genes$name
  colnames(counts) <- name(TF_clusters_PWMs$All_pwms_log_odds[sel])
  counts <- counts[rowSums(counts)>0,]
  counts <- t(counts)
  mygrid <- somgrid(xdim= 10, ydim= 10, topo = 'hexagonal', toroidal = T)
  set.seed(1)
  som.model <- supersom(counts, grid = mygrid, rlen = 100)
  saveRDS(som.model, "Rdata/som_informative_motifs.rds")
}
if(!file.exists("Rdata/motif_collapsed_cluster_counts.rds"))
{
  load("/groups/stark/almeida/data/motifs/motif_collection_v7_no_transfac_SteinAerts/TF_clusters_PWMs.RData")
  som <- readRDS("Rdata/som_informative_motifs.rds")
  dmat <- as.data.table(som$data[[1]], keep.rownames = T)
  colnames(dmat)[1] <- "motif"
  dmat[, cl:= som$unit.classif]
  mot <- melt(dmat, id.vars = c("cl", "motif"))
  mot[, Dmel:= .(TF_clusters_PWMs$metadata[match(motif, TF_clusters_PWMs$metadata$motif_name), "Dmel"])]
  mot[, Dmel:= paste(unique(na.omit(Dmel)), collapse = "_"), cl]
  mot <- mot[, .SD[motif==.SD[, sum(((2^value)-1)), motif][which.max(V1), motif]], cl]
  mot[, c("FBgn", "symbol", "gene_part"):= tstrsplit(variable, "_")]
  mot[grep("^TSS", gene_part), c("var", "value"):= .("TSS_count", log2(sum((2^value)-1)+1)), .(FBgn, symbol, motif)]
  mot[grep("^RE", gene_part), c("var", "value"):= .("RE_count", log2(sum((2^value)-1)+1)), .(FBgn, symbol, motif)]
  mot <- unique(mot[, .(motif, cl, value, FBgn, symbol, var, Dmel)])
  mot[, cl:= paste0("cl", cl)]
  res <- dcast(mot, FBgn+symbol~var+cl, value.var= "value", fill = 0)
  saveRDS(res, "Rdata/motif_collapsed_cluster_counts.rds")
  mot_dmel <- unique(mot[, .(motif, cl, Dmel)])
  saveRDS(mot_dmel, "Rdata/som_informative_motifs_Dmel_TF.rds")
}

#----------------------------------------------#
# 4- SOM FC
#----------------------------------------------#
if(!file.exists("Rdata/som_genes_FC_PH_vs_W.rds"))
{
  dat <- data.table(file= list.files("db/FC_tables_all_transcriptomes/", ".txt", full.names = T))
  dat[, exp := gsub("_FC.txt", "", basename(file))]
  dat <- dat[, fread(file), (dat)]
  dat <- dcast(dat, rn~exp, value.var = c("log2FoldChange", "padj"))
  dat <- dat[(abs(log2FoldChange_PH18_vs_W18)<1 & (padj_PH18_vs_W18>0.05 | is.na(padj_PH18_vs_W18)))]
  sub <- dat[, .(rn,
                 log2FoldChange_PHJ9_vs_WKD,
                 log2FoldChange_PHJ11_vs_WKD,
                 log2FoldChange_PH29_vs_W29)]
  mat <- as.matrix(sub, 1)
  size <- 28
  mygrid <- somgrid(xdim= size, ydim= size, topo = 'rectangular', toroidal = T)
  set.seed(1)
  som.model <- supersom(mat, grid = mygrid, rlen = 500)
  codes <- som.model$codes[[1]]
  codes <- apply(codes, 2, function(x)
  {
    lim <- quantile(x, c(0.01, 0.99))
    x[x<lim[1]] <- lim[1]
    x[x>lim[2]] <- lim[2]
    return(x)
  })
  set.seed(1)
  som.model$kcl <- kmeans(codes, centers = 8, algorithm = "Lloyd")$cluster
  som.model$unit.classif_kcl <- som.model$kcl[som.model$unit.classif]
  saveRDS(som.model, "Rdata/som_genes_FC_PH_vs_W.rds")
}

#----------------------------------------------#
# 4- Merge all data 
#----------------------------------------------#
# FC
FC <- data.table(file= list.files("db/FC_tables_all_transcriptomes/", ".txt", full.names = T))
FC[, exp := gsub("_FC.txt", "", basename(file))]
FC <- FC[, fread(file), (FC)]
FC <- dcast(FC, rn~exp, value.var = c("log2FoldChange", "padj"))
sel <- c("rn", "PCXT1092A_vs_2A", "PSCSUZ21B842D_vs_42D", "EZ7312A_vs_2A", "SUZ1212A_vs_2A",
         "PH18_vs_PH29", "PHJ9_vs_PH29", "PHJ11_vs_PH29",
         "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18",
         "PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29", "PHRNAI_vs_WRNAI",
         "WTE1416_vs_120hED", "72hED_vs_120hED", "96hED_vs_120hED",
         "72hED_vs_WTE1416", "96hED_vs_72hED", "120hED_vs_96hED")
FC <- FC[, grep(paste(sel, collapse= "|"), colnames(FC)), with= F]

# SOM clustering
som <- readRDS("Rdata/som_genes_FC_PH_vs_W.rds")
cl <- data.table(FBgn= rownames(som$data[[1]]))
cl[, c("som_cl", "som_kcl"):= .(som$unit.classif, som$unit.classif_kcl)]
dat <- merge(cl, FC, by.x= "FBgn", by.y= "rn")

# ChIP
ChIP <- readRDS("Rdata/ChIP_normalized_enrichment_TSS_REs.rds") 
dat <- merge(dat, ChIP, by= "FBgn", all.x= T)

# Motifs
mot <- readRDS("Rdata/motif_collapsed_cluster_counts.rds")
colnames(mot) <- gsub("_cl", "_motif_cl", colnames(mot))
dat <- merge(dat, mot, by= c("FBgn", "symbol"), all.x= T)

### SAVE
saveRDS(dat, "Rdata/genes_features_final.rds")


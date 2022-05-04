setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_RNA_CUTNRUN"]
meta <- meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), .(DESeq2_object, cdition)]
meta <- unique(meta[FC_file!="NA"])
meta <- meta[, fread(FC_file), (meta)]
sel <- meta[, any(padj<0.05 & abs(log2FoldChange)>1), FBgn][(V1), FBgn] # any signif foldChange
sel <- sel[!sel %in% meta[cdition=="RNA_PH18" & padj<0.05, FBgn]] # Not affected in ph18
meta <- meta[sel, on= "FBgn"][cdition!="RNA_PH18"]
meta <- meta[, cdition:= factor(cdition, c("RNA_PHD11", "RNA_PHD9", "RNA_PH29"))]
setorderv(meta, "cdition")

# Clip outliers
meta[, corr_log2FoldChange:= {
  lim <- quantile(log2FoldChange, c(0.05, 0.95), na.rm= T)
  log2FoldChange[log2FoldChange<lim[1]] <- lim[1]
  log2FoldChange[log2FoldChange>lim[2]] <- lim[2]
  scale(log2FoldChange)
}, FC_file]
meta[padj>0.05, c("log2FoldChange", 
                  "corr_log2FoldChange"):= .(NA, 0)]

# Clustering
layers_NA <- as.matrix(dcast(meta, FBgn~cdition, value.var = "log2FoldChange"), 1)
layers <- as.matrix(dcast(meta, FBgn~cdition, value.var = "corr_log2FoldChange"), 1)
layers <- list(constant= layers[, "RNA_PH29", drop= F],
               transient= layers[, c("RNA_PHD11", "RNA_PHD9")])
grid <- somgrid(2, 
                3, 
                "hexagonal", 
                toroidal= T)
init <- lapply(layers, function(x)
{
  set.seed(8)
  x <- x[sample(nrow(x), grid$xdim*grid$ydim), , drop= F]
  return(x)
})
som <- supersom(data = layers, 
                grid= grid,
                init = init,
                user.weights= c(1,1), 
                maxNA.fraction = 1)
som$ordered_matrices <- layers_NA[order(som$unit.classif),]

# Cluster object
cl <- list(log2FC_mat= layers_NA,
           corrLog2FC_mat= layers,
           som= som,
           rows= data.table(FBgn= rownames(som$data[[1]]), 
                            order= order(som$unit.classif), 
                            cl= som$unit.classif))

# Add gene names
FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
FBGN <- as.data.table(FBGN)
cl$rows[FBGN, symbol:= i.gene_symbol, on= "FBgn==gene_id"]
cl$rows[, col:= vl_palette_few_categ(max(cl))[cl]]

# Add PRC1 binding
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
PRC1 <- loadRData("external_data/SA2020_cl.list")
cl$rows[, PRC1_bound:= symbol %in% unlist(PRC1$genes)]

# Add FPKMs
fpkms <- as.data.table(DESeq2::fpkm(readRDS("db/dds/RNA/epiCancer_ED_RNA_CUTNRUN_dds.rds")), 
                       keep.rownames = "FBgn")
fpkms <- melt(fpkms, id.vars = "FBgn")
fpkms[, variable:= gsub("_1.bam|_2.bam|_3.bam", "_FPKM", variable)]
fpkms <- dcast(fpkms, FBgn~variable, value.var = "value", fun.aggregate = mean)
cl$rows <- cbind(cl$rows, fpkms[cl$rows$FBgn, on="FBgn"])

# Add HTMs
regions <- FBGN[type=="gene"][cl$rows, .(seqnames, start, end, cl), on= "gene_id==FBgn"]
regions <- vl_resizeBed(regions, upstream = 2500, downstream = 2500, center = "center")
regions[, H3K27me3_W18:= vl_bw_coverage(regions, 
                                        "db/bw/cutnrun/H3K27me3_PH18_merge.bw")]
regions[, H3K27Ac_W18:= vl_bw_coverage(regions, 
                                       "db/bw/cutnrun/H3K27Ac_PH18_merge.bw")]
cl$HTMs <- regions

# Add chromatin types SA2020
chromhmm <- vl_importBed("external_data/chromatin_types_SA2020_table_s1.txt")
TSSs <- as.data.table(GRanges(cl$rows$TSS))
cl$rows$chromhmm <- chromhmm[TSSs, name, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$name


# Gene Ontologies
cl$GO <- vl_GO_clusters(FBgn_list = split(cl$rows$FBgn, cl$rows$cl), 
                        go_object = vl_Dmel_GO_FB2020_05,
                        plot= F)

# Gene network
cl$net <- vl_STRING_interaction(symbols = cl$rows$symbol, 
                                size = apply(do.call(cbind, cl$corrLog2FC_mat), 1, function(x) max(abs(x), na.rm= TRUE)),
                                plot= F, 
                                col = cl$rows$col)

# Motif enrichment at promoters
proms <- vl_resizeBed(FBGN[type=="gene"],
                      center = "start", 
                      upstream = 0, 
                      downstream = 0, 
                      ignore.strand = F, 
                      genome= "dm6")
ext_proms <- vl_resizeBed(proms,
                          center = "start", 
                          upstream = 500, 
                          downstream = 250, 
                          ignore.strand = F, 
                          genome= "dm6")
cl$rows[proms, TSS:= paste0(seqnames, ":", start, "-", end, ":", strand), on= "FBgn==gene_id"]
cl$rows[ext_proms, c("seqnames", "start", "end", "strand"):= .(seqnames, start, end, strand), on= "FBgn==gene_id"]
cl$rows[, seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm6, 
                                              names= seqnames, 
                                              start= start, 
                                              end= end, 
                                              strand= strand))]
mot_counts <- vl_motif_counts(cl$rows$seq)
colnames(mot_counts) <- paste0(colnames(mot_counts), "_mot")
cl$rows <- cbind(cl$rows,
                 mot_counts)
motif_cols <- colnames(mot_counts)
cl$motif_enr <- vl_motif_cl_enrich(as.matrix(cl$rows[, c("FBgn", motif_cols), with= F], 1), 
                                   cl_IDs = cl$rows$cl,
                                   plot= F)

# SAVE
saveRDS(cl,
        "Rdata/clustering_cutnrun_genotype_transcriptomes.rds")

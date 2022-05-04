setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA"]
meta <- meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), .(DESeq2_object, cdition)]
meta <- unique(meta[FC_file!="NA" & !grepl("RNA_PHD11.txt$", FC_file)]) # remove NAs and PHD11T vs PHD11
meta <- meta[cdition %in% c("RNA_PH18", "RNA_PHD11_T2", "RNA_PHD11_T3")]
meta <- meta[, fread(FC_file), (meta)]
sel <- meta[, any(padj<0.05 & abs(log2FoldChange)>1), FBgn][(V1), FBgn] # any signif foldChange
sel <- sel[!sel %in% meta[cdition=="RNA_PH18" & padj<0.05, FBgn]] # Not affected in ph18
dat <- meta[sel, on= "FBgn"][cdition != "RNA_PH18"]
dat <- dat[, cdition:= factor(cdition, c("RNA_PHD11_T2", "RNA_PHD11_T3"))]
setorderv(dat, "cdition")

# Clip outliers
dat[, corr_log2FoldChange:= {
  lim <- quantile(log2FoldChange, c(0.05, 0.95), na.rm= T)
  log2FoldChange[log2FoldChange<lim[1]] <- lim[1]
  log2FoldChange[log2FoldChange>lim[2]] <- lim[2]
  scale(log2FoldChange)
}, FC_file]
dat[padj>0.05, c("log2FoldChange", 
                 "corr_log2FoldChange"):= .(NA, 0)]

# Clustering
layers_NA <- as.matrix(dcast(dat, FBgn~cdition, value.var = "log2FoldChange"), 1)
layers <- as.matrix(dcast(dat, FBgn~cdition, value.var = "corr_log2FoldChange"), 1)
grid <- somgrid(2, 
                3, 
                "hexagonal", 
                toroidal= T)
set.seed(1)
init <- layers[sample(nrow(layers), grid$xdim*grid$ydim), , drop= F]
som <- supersom(data = layers, 
                grid= grid,
                init = init,
                maxNA.fraction = 1)

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

# Gene Ontologies
cl$GO <- vl_GO_clusters(FBgn_list = split(cl$rows$FBgn, cl$rows$cl), 
                        go_object = vl_Dmel_GO_FB2020_05,
                        plot= F)

# Gene network
cl$net <- vl_STRING_interaction(symbols = cl$rows$symbol, 
                                size = apply(cl$corrLog2FC_mat, 1, function(x) max(abs(x), na.rm= TRUE)),
                                plot= F, 
                                col = cl$rows$col)

# Motif enrichment at promoters
cl$motif_enr$seq <- vl_resizeBed(FBGN[type=="gene"][cl$rows, .(seqnames, start, end, strand), on= "gene_id==FBgn"],
                                 center = "start", 
                                 upstream = 500, 
                                 downstream = 250, 
                                 ignore.strand = F, 
                                 genome= "dm6")
cl$motif_enr$seq[, seq:= as.character(
  BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm6, 
                   names= seqnames, 
                   start= start, 
                   end= end, 
                   strand= strand))]

cl$motif_enr$counts <- vl_motif_counts(cl$motif_enr$seq$seq)

cl$motif_enr$enr <- vl_motif_cl_enrich(cl$motif_enr$counts, 
                                       cl_IDs = cl$rows$cl,
                                       plot= F)

#FPKMs
fpkms <- as.data.table(DESeq2::fpkm(readRDS("db/dds/RNA/epiCancer_ED_allograft_RNA_gDNA_dds.rds")), 
                       keep.rownames = "FBgn")
cl$ref_fpkms <- fpkms[, .(FBgn,
                          PH29= rowMeans(.SD[, .(TPH29_1.bam, TPH29_2.bam, TPH29_3.bam)]))]

#PRC1 binding
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
PRC1 <- loadRData("external_data/SA2020_cl.list")
cl$rows[, PRC1_bound:= symbol %in% unlist(PRC1$genes)]

# HTMs
regions <- FBGN[type=="gene"][cl$rows, .(seqnames, start, end, cl), on= "gene_id==FBgn"]
regions <- vl_resizeBed(regions, 5000, center = "center")
regions[, H3K27me3_W18:= vl_bw_coverage(regions, 
                                        "db/bw/cutnrun/H3K27me3_PH18_merge.bw")]
regions[, H3K27Ac_W18:= vl_bw_coverage(regions, 
                                       "db/bw/cutnrun/H3K27Ac_PH18_merge.bw")]
cl$HTMs <- regions

# SAVE
saveRDS(cl,
        "Rdata/clustering_allograft_transcriptomes.rds")











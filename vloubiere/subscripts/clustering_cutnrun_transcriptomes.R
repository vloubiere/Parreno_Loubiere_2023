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

# SAVE
saveRDS(cl,
        "Rdata/clustering_cutnrun_genotype_transcriptomes.rds")

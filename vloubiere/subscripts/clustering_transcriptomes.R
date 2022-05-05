setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)
require(readxl)
require(GenomicRanges)

obj <- list()
for(cdition in c("cutnrun", "allograft"))
{
  if(cdition=="cutnrun")
  {
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
  }else if(cdition=="allograft")
  {
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
  }
  # Cluster object
  cl <- list(log2FC_mat= layers_NA,
             corrLog2FC_mat= layers,
             som= som,
             rows= data.table(FBgn= rownames(som$data[[1]]), 
                              order= order(som$unit.classif), 
                              cl= som$unit.classif))
  
  # Add gene names and TSS coordinates
  FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
  GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
  FBGN <- as.data.table(FBGN)[type=="gene"]
  cl$rows[FBGN, symbol:= .(i.gene_symbol), on= "FBgn==gene_id"]
  cl$rows[FBGN, TSS:= paste0(i.seqnames, ":", ifelse(strand=="-", end, start), ":", strand), on= "FBgn==gene_id"]
  cl$rows[, promoter:= vl_resizeBed(GRanges(TSS),
                                    upstream = 500, 
                                    downstream = 250, 
                                    ignore.strand = F)[, paste0(seqnames, ":", start, "-", end, ":", strand)]]
  cl$rows[, ext_prom:= vl_resizeBed(GRanges(TSS),
                                    upstream = 2500, 
                                    downstream = 2500, 
                                    ignore.strand = F)[, paste0(seqnames, ":", start, "-", end, ":", strand)]]
  cl$rows[, col:= vl_palette_few_categ(max(cl))[cl]]
  
  # Add PRC1 binding
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  PRC1 <- loadRData("external_data/SA2020_cl.list")
  PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
  cl$rows[PRC1, PRC1_cluster:= i.PRC1_cluster, on= "symbol"]
  
  # Add FPKMs
  fpkms <- as.data.table(DESeq2::fpkm(readRDS("db/dds/RNA/epiCancer_ED_RNA_CUTNRUN_dds.rds")), 
                         keep.rownames = "FBgn")
  fpkms <- melt(fpkms, id.vars = "FBgn")
  fpkms[, variable:= gsub("_1.bam|_2.bam|_3.bam", "_FPKM", variable)]
  fpkms <- dcast(fpkms, FBgn~variable, value.var = "value", fun.aggregate = mean)
  cols <- names(fpkms)[-1]
  cl$rows[, (cols):= fpkms[.BY, cols, with=F], FBgn]
  
  # Add HTMs
  cl$rows[, H3K27me3_W18:= vl_bw_coverage(GRanges(ext_prom),
                                          "db/bw/cutnrun/H3K27me3_PH18_merge.bw")]
  cl$rows[, H3K27Ac_W18:= vl_bw_coverage(GRanges(ext_prom),
                                         "db/bw/cutnrun/H3K27Ac_PH18_merge.bw")]

  # Add chromatin types SA2020
  chromhmm <- vl_importBed("external_data/chromatin_types_SA2020_table_s1.txt")
  cl$rows$chromhmm <- chromhmm[as.data.table(GRanges(cl$rows$TSS)), name, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$name

  # Gene Ontologies
  cl$GO <- vl_GO_clusters(FBgn_list = split(cl$rows$FBgn, cl$rows$cl),
                          go_object = vl_Dmel_GO_FB2020_05,
                          plot= F)

  # Gene network
  if(cdition=="cutnrun")
    size <- apply(do.call(cbind, cl$corrLog2FC_mat), 1, function(x) max(abs(x), na.rm= TRUE))
    else if(cdition=="allograft")
      size <- apply(cl$corrLog2FC_mat, 1, function(x) max(abs(x), na.rm= TRUE))
  cl$net <- vl_STRING_interaction(symbols = cl$rows$symbol,
                                  size = size,
                                  plot= F,
                                  col = cl$rows$col)

  # Motif enrichment at promoters
  cl$rows[, seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm6,
                                                names= GRanges(promoter)))]
  mot_counts <- vl_motif_counts(cl$rows$seq)
  colnames(mot_counts) <- paste0(colnames(mot_counts), "_mot")
  cl$rows <- cbind(cl$rows,
                   mot_counts)
  motif_cols <- colnames(mot_counts)
  cl$motif_enr <- vl_motif_cl_enrich(as.matrix(cl$rows[, c("FBgn", motif_cols), with= F], 1),
                                     cl_IDs = cl$rows$cl,
                                     plot= F)
  obj[[cdition]] <- cl
}

# SAVE
saveRDS(obj,
        "Rdata/clustering_RNA.rds")

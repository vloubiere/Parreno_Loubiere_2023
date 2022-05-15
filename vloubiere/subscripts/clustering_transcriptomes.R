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
  dat <- merge(as.data.table(layers_NA, keep.rownames = "FBgn"),
               data.table(FBgn= rownames(som$data[[1]]), 
                          order= order(som$unit.classif), 
                          cl= som$unit.classif),
                         by= "FBgn")
  
  # Add gene names and TSS coordinates
  FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
  GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
  FBGN <- as.data.table(FBGN)[type=="gene"]
  dat[FBGN, symbol:= .(i.gene_symbol), on= "FBgn==gene_id"]
  dat[FBGN, TSS:= paste0(i.seqnames, ":", ifelse(strand=="-", end, start), ":", strand), on= "FBgn==gene_id"]
  dat[, promoter:= vl_resizeBed(GRanges(TSS),
                                    upstream = 500, 
                                    downstream = 250, 
                                    ignore.strand = F)[, paste0(seqnames, ":", start, "-", end, ":", strand)]]
  dat[, ext_prom:= vl_resizeBed(GRanges(TSS),
                                    upstream = 2500, 
                                    downstream = 2500, 
                                    ignore.strand = F)[, paste0(seqnames, ":", start, "-", end, ":", strand)]]
  dat[, col:= vl_palette_few_categ(max(cl))[cl]]
  
  # Add PRC1 binding
  loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
  }
  PRC1 <- loadRData("external_data/SA2020_cl.list")
  PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
  dat[PRC1, PRC1_cluster:= i.PRC1_cluster, on= "symbol"]
  
  # Add FPKMs
  fpkms <- as.data.table(DESeq2::fpkm(readRDS("db/dds/RNA/epiCancer_ED_RNA_CUTNRUN_dds.rds")), 
                         keep.rownames = "FBgn")
  fpkms <- melt(fpkms, id.vars = "FBgn")
  fpkms[, variable:= gsub("_1.bam|_2.bam|_3.bam", "_FPKM", variable)]
  fpkms <- dcast(fpkms, FBgn~variable, value.var = "value", fun.aggregate = mean)
  cols <- names(fpkms)[-1]
  dat[, (cols):= fpkms[.BY, cols, with=F], FBgn]
  
  # Add HTMs
  files <- list.files("db/bw/cutnrun/", "H3K27.*_merge.bw", full.names = T)
  cols <- gsub("_merge.bw", "", basename(files))
  dat[, (cols):= lapply(files, function(x) vl_bw_coverage(GRanges(ext_prom), x))]
  
  # Add chromatin types SA2020
  chromhmm <- vl_importBed("external_data/chromatin_types_SA2020_table_s1.txt")
  dat$chromhmm <- chromhmm[as.data.table(GRanges(dat$TSS)), name, .EACHI, on= c("seqnames", "start<=end", "end>=start")]$name

  # Gene Ontologies
  GO <- vl_GO_clusters(FBgn_list = split(dat$FBgn, dat$cl),
                       go_object = vl_Dmel_GO_FB2020_05,
                       plot= F)

  # Gene network
  if(cdition=="cutnrun")
    size <- apply(do.call(cbind, layers), 1, function(x) max(abs(x), na.rm= TRUE))
    else if(cdition=="allograft")
      size <- apply(layers, 1, function(x) max(abs(x), na.rm= TRUE))
  net <- vl_STRING_interaction(symbols = dat$symbol,
                               size = size,
                               plot= F,
                               col = dat$col)

  # Motif enrichment at promoters
  dat[, seq:= as.character(BSgenome::getSeq(BSgenome.Dmelanogaster.UCSC.dm6,
                                            names= GRanges(promoter)))]
  mot_counts <- vl_motif_counts(dat$seq)
  colnames(mot_counts) <- paste0(colnames(mot_counts), "_mot")
  rownames(mot_counts) <- dat$FBgn
  motif_enr <- vl_motif_cl_enrich(mot_counts,
                                  cl_IDs = dat$cl,
                                  plot= F)
  cl <- list(data= dat,
             motif_counts= mot_counts,
             mot_enrichment= motif_enr,
             GO= GO,
             network= net,
             som= som)
  obj[[cdition]] <- cl
}

# SAVE
saveRDS(obj,
        "Rdata/clustering_RNA.rds")

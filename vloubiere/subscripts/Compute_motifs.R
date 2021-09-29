RE <- readRDS("Rdata/gene_REs.rds")
check <- seqlengths(BSgenome.Dmelanogaster.UCSC.dm6)
check <- data.table(seqnames= names(check), 
                    length= check)
RE_sel <- as.data.table(GRanges(RE$RE_coor))
RE_sel[check, check := end<=i.length, on= c("seqnames")]
RE <- RE[(RE_sel$check)]

# Count motifs
if(!file.exists("Rdata/all_motif_counts.rds"))
{
  hit <- matchMotifs(vl_Dmel_motifs_DB$All_pwms_log_odds, 
                     GRanges(RE$RE_coor), 
                     genome= "dm6", 
                     p.cutoff= 5e-4, 
                     bg= "even", 
                     out= "scores")
  counts <- as.matrix(motifCounts(hit))
  colnames(counts) <- name(vl_Dmel_motifs_DB$All_pwms_log_odds)
  counts <- as.data.table(counts)
  res <- cbind(RE, counts)
  
  # Add clusters
  cl <- readRDS("Rdata/final_clustering_transcriptomes.rds")
  cl <- unique(cl[, .(FBgn, cl)])
  res[cl, cl:= i.cl, on= "FBgn"]
  res[is.na(cl), cl:= 0]
  setcolorder(res, "cl")
  
  saveRDS(res, "Rdata/all_motif_counts.rds")
}

if(!file.exists("Rdata/final_motifs_cl_enrich.rds"))
{
  res <- readRDS("Rdata/all_motif_counts.rds")
  cols <- colnames(res)[-c(1:7)]
  
  # All REs
  enrich <- list()
  for(i in seq(max(res$cl)))
  {
    .c <- lapply(cols, function(x)
    {
      fisher.test(table(res$cl==i, res[[x]]>0))
    })
    
    padj <- p.adjust(sapply(.c, function(x) x$p.value), method= "fdr")
    log2FE <- log2(sapply(.c, function(x) x$estimate))
    enrich[[i]] <- data.table(cl= i,
                              motif= cols,
                              padj= padj,
                              log2FE= log2FE)
  }
  all <- rbindlist(enrich)
  
  # CPs
  enrich <- list()
  res <- readRDS("Rdata/all_motif_counts.rds")
  res <- res[grepl("^TSS", RE_ID)]
  for(i in seq(max(res$cl)))
  {
    .c <- lapply(cols, function(x)
    {
      fisher.test(table(res$cl==i, res[[x]]>0))
    })
    
    padj <- p.adjust(sapply(.c, function(x) x$p.value), method= "fdr")
    log2FE <- log2(sapply(.c, function(x) x$estimate))
    enrich[[i]] <- data.table(cl= i,
                              motif= cols,
                              padj= padj,
                              log2FE= log2FE)
  }
  CPs <- rbindlist(enrich)
  
  # Enhancers
  enrich <- list()
  res <- readRDS("Rdata/all_motif_counts.rds")
  res <- res[grepl("^RE", RE_ID)]
  for(i in seq(max(res$cl)))
  {
    .c <- lapply(cols, function(x)
    {
      fisher.test(table(res$cl==i, res[[x]]>0))
    })
    
    padj <- p.adjust(sapply(.c, function(x) x$p.value), method= "fdr")
    log2FE <- log2(sapply(.c, function(x) x$estimate))
    enrich[[i]] <- data.table(cl= i,
                              motif= cols,
                              padj= padj,
                              log2FE= log2FE)
  }
  ENs <- rbindlist(enrich)
  
  final <- list(all_REs= all, 
                TSSs= CPs, 
                enhancers= ENs)
  saveRDS(final, "Rdata/final_motifs_cl_enrich.rds")
}

# Selct motifs significantly enriched in each catergory
mot <- readRDS("Rdata/final_motifs_cl_enrich.rds")
meta <- as.data.table(vl_Dmel_motifs_DB$metadata)
mot <- lapply(mot, function(x)
{
  x[meta, symbol:= i.Dmel, on= "motif==motif_id"]
  x[, symbol2:= tstrsplit(symbol, "/", keep= 1), symbol] # Solves EcR/usp
  x[, sel:= any(padj<0.001)
      & all(is.finite(log2FE))
      & any(log2FE>1) 
      & !is.na(symbol), motif]
  x <- x[(sel)]
  symbols <- unique(import("../../genomes/dm6/dmel-all-r6.36.gtf")$gene_symbol)
  x[, sel := any(grepl(symbol, symbols)) | any(grepl(symbol2, symbols))
      , .(symbol, symbol2)]
  x <- x[(sel)]
  x <- x[, {
    .SD[motif==motif[which.min(padj)]]
    
  }, symbol]
  return(x[(sel), !"sel"])
})

saveRDS(mot, "Rdata/selected_motifs.rds")

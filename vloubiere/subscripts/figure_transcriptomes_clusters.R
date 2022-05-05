setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import data
obj <- readRDS("Rdata/clustering_RNA.rds")

pdf("pdf/Figures/Clustering_allograft_cutnrun_systems_RNA.pdf", 
    width= 25)
mat <- matrix(c(1,2,3,6,7,8,
                1,2,4,6,7,8,
                1,2,5,6,7,8),
              nrow= 3, byrow = T)
layout(mat, 
       widths = c(1,1,0.9,1.5,2.4,1))
for(cdition in names(obj))
{
  list2env(obj[[cdition]], 
           envir = environment())
  if(cdition=="allograft")
  {
    lims <- c(-8, 8)
    adj.vertices <- 10
    top_network <- 300
    mot_padj_cutoff <- 0.00001
    top_mot_enr <- 10
    cex_ball_GO <- 0.75
  }
  if(cdition=="cutnrun")
  {
    lims <- c(-5, 5)
    adj.vertices <- 4
    top_network <- 500
    mot_padj_cutoff <- 0.01
    top_mot_enr <- 6
    cex_ball_GO <- 0.5
  }
  
  # Heatmap
  par(mar= c(7,10,2,10))
  vl_heatmap(log2FC_mat[rows$order,], 
             cluster_rows= F,
             cluster_cols= F, 
             breaks = seq(lims[1], lims[2], length.out= 100), 
             show_rownames = F, 
             col = vl_palette_blueWhiteRed(100), 
             legend_title = "Fold Change (log2)", 
             auto_margins = F)
  cl_counts <- table(som$unit.classif)
  cl_pos <- cumsum(rev(cl_counts))
  abline(h= cl_pos[-length(cl_pos)], 
         lwd= 2,
         xpd= F)
  text(par("usr")[1], 
       cl_pos-diff(c(1, cl_pos))/2,
       rev(paste0("Cluster ", names(cl_counts), "\n(", cl_counts, ")")),
       xpd= T,
       pos= 2)
  title(paste0(cdition, " (SOM)"))
  
  # Histone marks
  mat <- apply(rows[(order), .(H3K27Ac_W18, H3K27me3_W18)], 2, function(x) log2(x+0.001))
  vl_heatmap(scale(mat), 
             cluster_rows= F,
             show_rownames= F,
             breaks= c(-2,0,2), 
             auto_margins = F)
  abline(h= cl_pos[-length(cl_pos)], 
         lwd= 2,
         xpd= F)
  text(par("usr")[1], 
       cl_pos-diff(c(1, cl_pos))/2,
       rev(paste0("Cluster ", names(cl_counts), "\n(", cl_counts, ")")),
       xpd= T,
       pos= 2)
  
  # CHROMHMM
  all_chrom <- fread("external_data/chromatin_types_SA2020_table_s1.txt")
  chrom <- rbind(rows[, .(chromhmm, cl)],
                 all_chrom[, .(chromhmm= V4, cl="all")])
  chrom <- dcast(na.omit(chrom), 
                 chromhmm~cl, 
                 fun.aggregate = length)
  chrom <- as.matrix(chrom, 1)
  chrom <- apply(chrom, 2, function(x) x/sum(x)*100)
  chrom <- chrom[order(apply(chrom, 1, sum), decreasing = T),]
  Cc <- sapply(rownames(chrom), function(x) switch(x, 
                                                   "aTSS"= "tomato",
                                                   "aTTS"= "red",
                                                   "Enhancer"= "limegreen",
                                                   "PcG"= "cornflowerblue",
                                                   "Null"= "grey"))
  par(mar= c(2,5,2,6))
  barplot(chrom,
          col= Cc,
          las= 1,
          ylab= "% chromhmm class")
  legend(par("usr")[2],
         par("usr")[4],
         legend = rownames(chrom),
         fill= Cc,
         xpd= T,
         bty= "n")
  
  # PRC1 binding
  all_PRC1 <- loadRData("external_data/SA2020_cl.list")
  all_PRC1 <- rbindlist(lapply(all_PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
  PRC1_binding <- rbind(rows[, .(PRC1_cluster, cl)],
                        all_PRC1[, .(PRC1_cluster, cl= "all")])
  PRC1_binding[, total:= ifelse(cl=="all", 16514, .N), cl]
  PRC1_binding <- PRC1_binding[!is.na(PRC1_cluster)]
  PRC1_binding[, PRC1_cluster:= switch(PRC1_cluster,
                                       "aTSS A"= "aTSS", 
                                       "aTSS B"= "aTSS", 
                                       "Enhancer A"= "Enhancer", 
                                       "Enhancer B"= "Enhancer", 
                                       "Polycomb A"= "Polycomb", 
                                       "Polycomb B"= "Polycomb", 
                                       "Repressed"= "Repressed"), PRC1_cluster]
  PRC1_binding <- PRC1_binding[, .(perc= .N/total*100), .(cl, PRC1_cluster)]
  mat <- dcast(unique(PRC1_binding), 
               PRC1_cluster~cl, 
               value.var = "perc", 
               fill = 0)
  mat <- as.matrix(mat, 1)
  mat <- mat[c("aTSS", "Enhancer", "Polycomb", "Repressed"),]
  Cc <- c("tomato", "limegreen", "cornflowerblue", "grey")
  bar <- barplot(mat,
                 ylab= "% PRC1 bound genes",
                 las= 1,
                 col= Cc)
  legend(par("usr")[2],
         par("usr")[4],
         legend = rownames(mat),
         fill= Cc,
         xpd= T,
         bty= "n")
  
  # ref FPKMs
  if(cdition=="allograft")
    .f <- as.formula(PH29_FPKM~cl) else if(cdition=="cutnrun")
      .f <- as.formula(W18_FPKM~cl)
  vl_boxplot(.f,
             rows,
             las= 0,
             ylab= .f[[2]])
  
  # Network
  par(mar= c(0,0,1,0),
      las= 0)
  plot(net,
       size= net$vertices$size*adj.vertices,
       top_N = top_network)
  leg <- unique(rows[net$vertices, .(cl, col), on= "symbol==name"])[order(cl)]
  legend("topleft",
         fill= leg$col,
         legend= paste0("Cluster ", leg$cl),
         bty= "n")
  title("STRING interactions")
  
  # GOs
  par(las= 1,
      mar= c(2,35,1,8))
  plot(GO, 
       padj_cutoff= 0.001, 
       auto_margins= F, 
       N_top= 10, 
       cex.balloons= cex_ball_GO)
  title("GOs enrichment")
  
  # Motifs
  par(las= 1,
      mar= c(2,8,2,5))
  plot(motif_enr, 
       padj_cutoff= mot_padj_cutoff, 
       auto_margins= F, 
       N_top= top_mot_enr,
       cex.balloons= 1.5)
  title("Motifs enrichment")
}
dev.off()
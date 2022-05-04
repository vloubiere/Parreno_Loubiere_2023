setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import data
obj <- lapply(c("Rdata/clustering_allograft_transcriptomes.rds",
                "Rdata/clustering_cutnrun_genotype_transcriptomes.rds"), function(x) readRDS(x))
names(obj) <- c("allograft", "cuntrun")

pdf("pdf/Figures/Clustering_transcriptomes.pdf", 
    width= 25)
mat <- matrix(c(1,2,3,5,6,7,
                1,2,4,5,6,7),
              nrow= 2, byrow = T)
layout(mat, widths = c(1,1,0.9,1.5,2.4,1))
for(i in seq(obj))
{
  list2env(obj[[i]], 
           envir = environment())
  if(i==1)
  {
    lims <- c(-8, 8)
    adj.vertices <- 10
    top_network <- 300
    mot_padj_cutoff <- 0.00001
    top_mot_enr <- 10
    cex_ball_GO <- 0.75
  }
  if(i==2)
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
  title(paste0(c("Allograft", "CutNrun genotype")[i], " (SOM)"))
  
  # Histone marks
  mat <- apply(HTMs[rows$order, .(H3K27Ac_W18, H3K27me3_W18)], 2, function(x) log2(x+0.001))
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
  
  # PRC1 binding
  par(mar= c(5,5,4,2))
  pl <- table(rows$PRC1_bound, rows$cl)
  perc <- c(genome= 2069/16514*100, # SA202 PRC1 bound genes vs all mRNA/ncRNA genes
            apply(pl, 2, function(x) round(x[2]/sum(x)*100)))
  bar <- barplot(perc,
                 ylab= "% PRC1 bound genes",
                 las= 1)
  text(bar, 
       perc,
       c(16514, apply(pl, 2, sum)),
       xpd= T,
       pos= 3)
  
  # ref FPKMs
  cols <- names(ref_fpkms)[-1]
  ref_fpkms[rows, cl:= i.cl, on= "FBgn"]
  ref_fpkms[, cl:= as.character(cl)]
  .c <- rbind(na.omit(ref_fpkms),
              ref_fpkms[, cl:="genome"])
  .c[, cl:= factor(cl)]
  .c[, cl:= relevel(cl, "genome")]
  par(las= 2)
  vl_boxplot(as.formula(paste0(cols, "~cl")),
             .c,
             las= 0,
             ylab= paste("FPKMs", cols))
  
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
  plot(motif_enr$enr, 
       padj_cutoff= mot_padj_cutoff, 
       auto_margins= F, 
       N_top= top_mot_enr,
       cex.balloons= 1.5)
  title("Motifs enrichment")
}
dev.off()
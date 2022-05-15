setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import data
obj <- readRDS("Rdata/clustering_RNA.rds")

pdf("pdf/Figures/Cluster_PRC1_bound_unbound_RNA.pdf",
    width= 9,
    height= 6)
par(las= 1,
    mar= c(2,3,4,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4))
# for(cdition in names(obj))
for(cdition in "cutnrun")
{
  list2env(obj[[cdition]], 
           envir = environment())
  for(cols in colnames(log2FC_mat))
  {
    .c <- split(log2FC_mat[,cols], rows[, .(ifelse(!is.na(PRC1_cluster), "PRC1", "noPRC1"), paste0("cl_", cl))])
    
    vl_boxplot(.c,
               ylab= "log2FoldChange", 
               boxcol= c("lightgrey", "tomato"),
               compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
               xaxt= "n")
    abline(h= 0, lty= "11")
    axis(1, 
         at = seq(1.5, 11.5, 2), 
         labels = paste0("cl ", 1:6))
    title(main= cols, 
          line=3)
    legend(6,
           par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
           legend= c("noPRC1", "PRC1 bound"),
           fill= c("lightgrey", "tomato"),
           bty= "n",
           xpd=T)
  }
  .c <- split(rows$PH18_FPKM, rows[, .(ifelse(!is.na(PRC1_cluster), "PRC1", "noPRC1"), paste0("cl_", cl))])
  vl_boxplot(.c,
             ylab= "FPKM PH18 (log2)", 
             xaxt= "n",
             compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
             boxcol= c("lightgrey", "tomato"))
  title(main= "FPKMs PH18",
        line=3)
  axis(1, 
       at = seq(1.5, 11.5, 2), 
       labels = paste0("cl ", 1:6))
  legend(6,
         par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
         legend= c("noPRC1", "PRC1 bound"),
         fill= c("lightgrey", "tomato"),
         bty= "n",
         xpd=T)
  for(cols in unlist(lapply(c("H3K27me3_", "H3K27Ac_"), function(x) paste0(x, c("W18", "PHD11", "PHD9", "PH29")))))
  {
    .c <- split(rows[[cols]], rows[, .(ifelse(!is.na(PRC1_cluster), "PRC1", "noPRC1"), paste0("cl_", cl))])
    if(grepl("H3K27me3", cols))
      yl <- c(0,30) else if(grepl("H3K27Ac", cols))
        yl <- c(0,9)
    vl_boxplot(.c,
               ylab= "CUT&RUN enrichment", 
               boxcol= c("lightgrey", "tomato"),
               compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
               xaxt= "n",
               ylim= yl)
    axis(1, 
         at = seq(1.5, 11.5, 2), 
         labels = paste0("cl ", 1:6))
    title(main= gsub("W18", "PH18", cols), 
          line=3)
    legend(6,
           par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
           legend= c("noPRC1", "PRC1 bound"),
           fill= c("lightgrey", "tomato"),
           bty= "n",
           xpd=T)
  }
}
dev.off()

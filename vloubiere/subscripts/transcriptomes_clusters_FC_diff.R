setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import data
obj <- readRDS("Rdata/clustering_RNA.rds")
list2env(obj, 
         envir = environment())
cols <- c("RNA_PHD11", 
          "RNA_PHD9", 
          "RNA_PH29", 
          "PH18_FPKM",
          "H3K27Ac_PH18",
          "H3K27Ac_PHD11",
          "H3K27Ac_PHD9",
          "H3K27Ac_PH29",
          "H3K27me3_PH18",
          "H3K27me3_PHD11",
          "H3K27me3_PHD9",
          "H3K27me3_PH29",
          "PH_ED",
          "PC_ED",
          "SUZ12_ED",
          "H2AK118Ub_ED",
          "H3K4me1_ED",
          "H3K27me2_ED")

pdf("pdf/Figures/Cluster_PRC1_bound_unbound_RNA.pdf",
    width= 9,
    height= 6)
par(las= 1,
    mar= c(2,3.5,4,0.5),
    mgp= c(2.5, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4))
for(col in cols)
{
  .c <- split(data[[col]], data[, .(ifelse(!is.na(PRC1_cluster), "PRC1", "noPRC1"), paste0("cl_", cl))])
  yl <- fcase(grepl("^RNA_", col), "log2FoldChange",
              grepl("^H3K27", col), "CUT&RUN Enrichment",
              grepl("FPKM$", col), "FPKM",
              grepl("_ED$", col), "ChIP-Seq Enrichment")
  lims <- NULL
  if(grepl("^H3K27Ac_", col))
    lims <- c(0, 8)
  if(grepl("^H3K27me3_", col))
    lims <- c(0, 27)
  if(is.null(lims))
    vl_boxplot(.c,
               ylab= yl, 
               boxcol= c("lightgrey", "tomato"),
               compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
               xaxt= "n")
    else
      vl_boxplot(.c,
                 ylab= yl, 
                 boxcol= c("lightgrey", "tomato"),
                 compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12)),
                 ylim= lims,
                 xaxt= "n")
  abline(h= 0, lty= "11")
  axis(1, 
       at = seq(1.5, 11.5, 2), 
       labels = paste0("cl ", 1:6))
  title(main= col, 
        line=3)
  legend(6,
         par("usr")[4]+diff(grconvertY(c(0,3), "line", "user")), 
         legend= c("noPRC1", "PRC1 bound"),
         fill= c("lightgrey", "tomato"),
         bty= "n",
         xpd=T)
}
dev.off()

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
    width= 7,
    height= 4)
par(las= 2,
    mar= c(5,3,2,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(2,3))
for(cdition in names(obj))
{
  list2env(obj[[cdition]], 
           envir = environment())
  for(cols in colnames(log2FC_mat))
  {
    .c <- split(log2FC_mat[,cols], rows[, .(ifelse(!is.na(PRC1_cluster), "PRC1", "noPRC1"), paste0("cl_", cl))])
    vl_boxplot(.c,
               ylab= "log2FoldChange", 
               tilt.names = T)
    title(main= paste(cols, cdition))
  }
}
dev.off()

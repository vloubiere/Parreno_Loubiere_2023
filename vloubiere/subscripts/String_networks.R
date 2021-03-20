setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(data.table)
require(kohonen)
require(pheatmap)
require(ontologyIndex)
require(STRINGdb)

dat <- readRDS("Rdata/som_clustering_transcriptomes.rds")
string_db <- STRINGdb$new(version="11", species= 7227, score_threshold= 200, input_directory="") # Species ID can be found on STRING website, 11 is the current version of STRING
current <- as.data.table(string_db$map(dat, "symbol", removeUnmappedRows = TRUE)) # TRhe second argument corresponds to the colname within the data.table containing gene symbols

pdf("pdf/string_cluster_networks.pdf")
current[, {
  if(length(unique(STRING_id))<200){
    string_db$plot_network(STRING_id)
    mtext(paste("Cluster #", cl), side= 1)
  }
  print(paste0(cl, " DONE!"))
}, cl]
dev.off()

# example1_mapped_pval05 <- string_db$add_diff_exp_color(subset(example1_mapped, padj<0.05), logFcColStr="log2FoldChange")
# payload_id <- string_db$post_payload(example1_mapped_pval05$STRING_id, colors=example1_mapped_pval05$color )
# string_db$plot_network(hits, payload_id=payload_id)

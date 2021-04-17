dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")

# Species ID can be found on STRING website, 11 is the current version of STRING
string_db <- STRINGdb$new(version="11", species= 7227, score_threshold= 200, input_directory="") 

# TRhe second argument corresponds to the colname within the data.table containing gene symbols
current <- as.data.table(string_db$map(dat, "symbol", removeUnmappedRows = TRUE)) 

pdf("pdf/string_cluster_networks.pdf")
current[, {
  .c <- .SD[order(abs(log2FoldChange), decreasing = T)][1:200, STRING_id]
  string_db$plot_network(.c)
  mtext(paste("Cluster #", cl), side= 1)
  print(paste0(cl, " DONE!"))
}, cl]
dev.off()

# example1_mapped_pval05 <- string_db$add_diff_exp_color(subset(example1_mapped, padj<0.05), logFcColStr="log2FoldChange")
# payload_id <- string_db$post_payload(example1_mapped_pval05$STRING_id, colors=example1_mapped_pval05$color )
# string_db$plot_network(hits, payload_id=payload_id)

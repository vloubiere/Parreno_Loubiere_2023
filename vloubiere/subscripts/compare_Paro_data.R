cl <- readRDS("Rdata/final_clustering_transcriptomes.rds")
Paro <- data.table(file= list.files("db/FC_tables/", "Paro", full.names= T))
Paro <- Paro[, fread(file), file]
Paro[cl, cl:= i.cl, on= "V1==FBgn"]
Paro[, cdition:= gsub("RNA_Paro_2018_RNA_|_vs_RNA_WTED_FC.txt", "", basename(file))]

mat <- as.matrix(dcast(Paro[!is.na(cl)], cl~cdition, value.var = "log2FoldChange", fun.aggregate = mean), 1)
vl_heatmap(mat, 
           cluster_rows = F, 
           display_numbers = T, 
           breaks = c(-8,0,8))

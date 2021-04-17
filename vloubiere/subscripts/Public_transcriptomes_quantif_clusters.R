dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")
dat <- unique(dat[, .(FBgn, cl, symbol, ycoor)])
dat <- dat[, .(file= c("db/FC_tables/RNA_development_RNA_72hED_vs_RNA_WTE1416_FC.txt",
                       "db/FC_tables/RNA_development_RNA_96hED_vs_RNA_72hED_FC.txt",
                       "db/FC_tables/RNA_development_RNA_120hED_vs_RNA_96hED_FC.txt")), dat]
dat <- dat[, { 
  .c <- fread(file)
  cbind(.SD, .c[match(FBgn, V1), baseMean:padj])
}, file]
dat[, cdition := gsub("RNA_development_|_FC.txt|_RNA|RNA_", "", basename(file))]
dat[, cdition := gsub("_vs_", "/", cdition)]

saveRDS(dat, "Rdata/Transcriptomes_quantif_clusters.rds")

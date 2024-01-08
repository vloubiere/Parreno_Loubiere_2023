setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[, .(FBgn, symbol, cluster= cl, class, PcG_bound)]
fwrite(dat,
       "db/extended_data_table/Extended_Data_Table_2.txt",
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)

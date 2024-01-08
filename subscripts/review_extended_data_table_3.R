setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

dat <- readRDS("Rdata/final_ATAC_table.rds")
dat <- dat[, .(ID, cluster= cl)]
fwrite(dat,
       "db/extended_data_table/Extended_Data_Table_3.txt",
       col.names = T,
       row.names = F,
       sep= "\t",
       quote= F,
       na= NA)

RE <- readRDS("Rdata/gene_REs.rds")

# ChIP
dat <- data.table(file= list.files("../../projects/public_data/dm6/bw/", full.names = T))[-28] # REmove rep1 K9
dat <- dat[!grepl("INPUT", file) & grepl("_ED|L3", file)]
dat[, cdition:= gsub("_ED_merge.bw|_merge|.bw$", "", basename(file))]
ChIP <- dat[, cbind(RE, vl_bw_coverage(GRanges(RE$RE_coor), bw = file)), (dat)]
ChIP <- ChIP[, .(score= mean(score)), .(FBgn, cdition)]

# Developmental transcriptomes
dev <- data.table(file= c("db/FC_tables/RNA_development_RNA_72hED_vs_RNA_WTE1416_FC.txt",
                          "db/FC_tables/RNA_development_RNA_96hED_vs_RNA_72hED_FC.txt",
                          "db/FC_tables/RNA_development_RNA_120hED_vs_RNA_96hED_FC.txt"))
dev <- dev[, fread(file), file]
colnames(dev)[colnames(dev)=="V1"] <- "FBgn"
dev[, cdition:= gsub("RNA_development_RNA_|_FC.txt", "", basename(file))]


res <- rbindlist(list(ChIP= ChIP[, .(FBgn, cdition, score)], 
                      Transcriptome= dev[, .(FBgn, cdition, score= log2FoldChange)]), 
                 idcol = T)
saveRDS(res, "Rdata/available_data.rds")

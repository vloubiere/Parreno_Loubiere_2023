RE <- readRDS("Rdata/gene_REs.rds")

# ChIP
dat <- data.table(file= list.files("../../public_data/dm6/bw/", full.names = T))
dat <- dat[!grepl("INPUT", file)]
dat[, cdition:= gsub("_merge.bw", "", basename(file))]
ChIP <- dat[, cbind(RE, vl_bw_coverage(GRanges(RE$RE_coor), bw = file)), (dat)]

# Tanscriptomes
dat <- data.table(file= c("db/FC_tables/RNA_development_RNA_72hED_vs_RNA_WTE1416_FC.txt",
                          "db/FC_tables/RNA_development_RNA_96hED_vs_RNA_72hED_FC.txt",
                          "db/FC_tables/RNA_development_RNA_120hED_vs_RNA_96hED_FC.txt"))
dat <- dat[, fread(file), file]
colnames(dat)[colnames(dat)=="V1"] <- "FBgn"

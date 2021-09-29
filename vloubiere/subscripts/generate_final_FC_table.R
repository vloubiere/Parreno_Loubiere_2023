setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- data.table(pat= c("epiCancer_dose_ED.*PH18_vs_RNA_W18",
                         "epiCancer_dose_ED.*PH21_vs_RNA_W21",
                         "epiCancer_dose_ED.*PH25_vs_RNA_W25",
                         "epiCancer_dose_ED.*PH29_vs_RNA_W29",
                         "epiCancer_ED.*EZ18_vs_RNA_W18",
                         "epiCancer_ED.*EZD11_vs_RNA_WKD",
                         "epiCancer_ED.*EZD9_vs_RNA_WKD",
                         "epiCancer_ED.*EZ29_vs_RNA_W29",
                         "epiCancer_ED.*PH18_vs_RNA_W18",
                         "epiCancer_ED.*PHD11_vs_RNA_WKD",
                         "epiCancer_ED.*PHD9_vs_RNA_WKD",
                         "epiCancer_ED.*PH29_vs_RNA_W29"))
dat[, file:= list.files("db/FC_tables/", 
                        pat, 
                        full.names = T), pat]
dat[, cdition:= gsub(".*_RNA_(.*)_vs_RNA_(.*).txt$", "\\1", file), file]
dat[, cdition:= paste0(cdition, ifelse(grepl("epiCancer_dose", file), "_dose", "_TS")), cdition]
dat[, cdition:= factor(cdition, levels = dat$cdition)]
dat <- dat[, fread(file), (dat)]
colnames(dat)[4] <- "FBgn"
symbols <- fread("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36_gene_symbols.txt")
dat[symbols, symbol:= i.gene_symbol, on= "FBgn==gene_id"]
setcolorder(dat, c("file", "cdition", "FBgn", "symbol"))
dat <- dat[, !c("file", "pat", "lfcSE")]
saveRDS(dat, "Rdata/final_FC_table.rds")
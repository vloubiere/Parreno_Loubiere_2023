setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import K27me3 domains ----
K27 <- vl_importBed("db/bed/merged_K27_domains/PH18.bed")

# Import genes table ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat[!(class %in% c("Reversible", "Irreversible")), class:= "Other"]

# Define K27me3 domains containing reversible/irreversible genes ----
K27[, irreversible:= vl_covBed(K27, dat[class=="Irreversible"])]
K27[, reversible:= vl_covBed(K27, dat[class=="Reversible"])]
K27[, other:= vl_covBed(K27, dat[class=="Other"])]
K27[, irrev_symbol:= vl_intersectBed(dat[class=="Irreversible"], .SD)[, paste0(sort(symbol), collapse= ",")], seq(nrow(K27))]
K27[, rev_symbol:= vl_intersectBed(dat[class=="Reversible"], .SD)[, paste0(sort(symbol), collapse= ",")], seq(nrow(K27))]
K27[, other_symbol:= vl_intersectBed(dat[class=="Other"], .SD)[, paste0(sort(symbol), collapse= ",")], seq(nrow(K27))]
cols <- c("irrev_symbol", "rev_symbol", "other_symbol")
K27[, (cols):= lapply(.SD, function(x) ifelse(x=="", NA, x)), .SDcols= cols]
K27 <- K27[order(-irreversible, -(reversible-irreversible))]

# Save
saveRDS(K27,
        "Rdata/K27_domains_classif_reversible_irreversible.rds")
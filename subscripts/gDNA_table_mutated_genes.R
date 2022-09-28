setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!(PH18) & !is.na(type)]
dat$FBtr <- NULL
res <- dat[, lapply(.SD, unlist), setdiff(names(dat), c("FBgn", "symbol")), .SDcols= c("FBgn", "symbol")]
res <- res[type %in% c("nonsynonymous SNV", "stopgain")
           & class=="Tumor specific"]
counts <- dcast(res[, .(FBgn, symbol, alt_class, cdition)],
                FBgn+symbol~alt_class+cdition, 
                fun.aggregate = length)
counts <- na.omit(counts)
counts[, total_mutations:= rowSums(.SD), .SDcols= setdiff(names(counts), c("FBgn", "symbol"))]
GO <- readRDS("Rdata/GO_mutated_genes.rds")
GO <- GO[variable %in% c("actin cytoskeleton organization",
                         "contractile fiber",
                         "extracellular matrix structural constituent")]
GO[, variable:= gsub(" ", "_", variable), variable]
GO[, variable:= gsub("_organization|_structural_constituent", "", variable), variable]
GO[, variable:= gsub(" ", "_", variable), variable]
GO[, {
  col <- as.character(variable)
  .c <- FBgn
  counts[, (col):= FBgn %in% .c]
  print("")
}, variable]
setorderv(counts, "total_mutations", -1)
setcolorder(counts, 
            c("FBgn",
              "symbol", 
              "total_mutations",
              "actin_cytoskeleton",
              "contractile_fiber",
              "extracellular_matrix"))

fwrite(counts, 
       "Rdata/mutated_gene_counts.txt", 
       sep= "\t", 
       col.names = T, 
       na= NA)

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# pdf and data file
dat <- data.table(file= c(list.files("db", recursive = T, full.names = T),
                          list.files("pdf", recursive = T, full.names = T)))
dat[, dir:= dirname(file)]
dat[dir %in% c("db/peaks/ATAC", "db/peaks/cutnrun"), file:= ifelse(grepl("confident", file), file, NA)]
dat[dir %in% c("db/bed/ATAC/", "db/bed/cutnrun_EcR/"), file:= NA]
dat <- na.omit(dat)
dat[, dir.create(paste0("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/", dir), recursive = T), dir]
dat[, file.copy(normalizePath(file), 
                paste0("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/", file),
                overwrite = T)]

# master tables
file.copy("Rdata/final_gene_features_table.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/db/final_gene_features_table.txt", 
          overwrite = T)
file.copy("Rdata/list_genes_interest_heatmaps.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/db/list_genes_interest_heatmaps.txt",
          overwrite = T)
file.copy("Rdata/mutated_gene_counts.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/db/mutated_gene_counts.txt",
          overwrite = T)
file.copy("Rdata/genes_overlapping_highly_mutated_regions.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/db/genes_overlapping_highly_mutated_regions.txt",
          overwrite = T)
file.copy("git_epiCancer/epiCancer_presentation.html",
          paste0("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/Presentations/", Sys.Date(), "_presentation.html"),
          overwrite = T)
file.copy("git_epiCancer/styles.css",
          paste0("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/Presentations/styles.css"),
          overwrite = T)

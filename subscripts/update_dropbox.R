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
dat[, dir.create(paste0("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/", dir), 
                 recursive = T, 
                 showWarnings = F), dir]
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

# Tables for collaborators
dir.create("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/tables_collaborators/", showWarnings = F)
dat <- fread("Rdata/final_gene_features_table.txt")
cols <- c("FBgn", "symbol", 
          "log2FoldChange_PH18", "padj_PH18",           
          "log2FoldChange_PHD11", "padj_PHD11", 
          "log2FoldChange_PHD9", "padj_PHD9",
          "log2FoldChange_PH29", "padj_PH29")
fwrite(dat[recovery=="Recovery", ..cols],
       "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/tables_collaborators/Recovery_genes.txt", 
       col.names = T, 
       sep= "\t",
       na= NA,
       quote= F)
fwrite(dat[recovery=="noRecovery", ..cols],
       "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/tables_collaborators/noRecovery_genes.txt", 
       col.names = T, 
       sep= "\t",
       na= NA,
       quote= F)
fwrite(dat[cl==2, ..cols],
       "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/tables_collaborators/nonRecovering_cluster2_genes.txt", 
       col.names = T, 
       sep= "\t",
       na= NA,
       quote= F)
fwrite(dat[cl==5, ..cols],
       "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/tables_collaborators/Recovering_cluster5_genes.txt", 
       col.names = T, 
       sep= "\t",
       na= NA,
       quote= F)







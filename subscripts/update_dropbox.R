setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# 1st colum -> local folder | 2nd -> remote folder
dat <- matrix(c("db/dds/RNA/", # dds RNA
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/RNA/dds/",
                "db/dds/cutnrun/", # dds CUTNRUN
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/CUTNRUN/dds/",
                "db/FC_tables/RNA/", # FC RNA
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/RNA/FC_tables/",
                "db/FC_tables/cutnrun/", # FC CUTNRUN
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/CUTNRUN/FC_tables/",
                "db/bw/RNA-Seq/", # bw RNA
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/RNA/bw/",
                "db/bw/cutnrun/", # bw CUTNRUN
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/CUTNRUN/bw/",
                "db/peaks/cutnrun_merged_peaks/", # merged peaks CUTNRUN
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/CUTNRUN/merged_peaks_changes/", 
                "db/peaks/cutnrun/", #  peaks CUTNRUN
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/CUTNRUN/peaks/", 
                "pdf/", # All PDFs
                "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/pdf/"), 
              ncol= 2, 
              byrow = T)
dat <- as.data.table(dat)
setnames(dat, c("local", "dropbox"))
dat[, {
  file.remove(list.files(dropbox, 
                         full.names = T, 
                         recursive= T))
  file.copy(list.files(local, 
                       full.names = T, 
                       recursive= T),
            dropbox, 
            recursive = T)
  print(".")
}, (dat)]

# Indidivual files
file.copy("db/FC_tables/RNA_table_AMM.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/RNA/RNA_table_AMM.txt", overwrite = T)
file.copy("Rdata/RECOVERY_NORECOVERY_genes.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/RNA/RECOVERY_NORECOVERY_genes.txt", overwrite = T)
file.copy("db/FC_tables/CUTNRUN_table_AMM.txt",
          "/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/CUTNRUN/CUTNRUN_table_AMM.txt", overwrite = T)

file.copy("git_epiCancer/presentation.html",
          paste0("/mnt/c/Users/User/Dropbox (Compte personnel)/collaborations/epigenetic_cancer/Presentations/presentation_", Sys.Date(), ".html"), overwrite = T)

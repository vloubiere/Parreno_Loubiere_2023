setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)
require(readxl)

meta <- as.data.table(read_xlsx("Rdata/metadata_DNA.xlsx"))
comp <- data.table(comparison= c("PHD9_vs_PH18",
                                 "PHD11_vs_PH18",
                                 "PH29_vs_PH18",
                                 "HOST_vs_PH18",
                                 "PH29_T1_vs_primary",
                                 "PH29_T5_vs_primary",
                                 "PH29_T20_vs_primary",
                                 "PH29_T25_vs_primary",
                                 "PHD11_T20_vs_primary",
                                 "PHD11_T25_vs_primary",
                                 "PH29_T1_vs_host",
                                 "PH29_T5_vs_host",
                                 "PH29_T20_vs_host",
                                 "PH29_T25_vs_host",
                                 "PHD11_T20_vs_host",
                                 "PHD11_T25_vs_host"), 
                   file_pattern=c("PH18_1_PHD9_1.somaticcall",
                                  "PH18_1_PHD11_1.somaticcall",
                                  "PH18_1_PH29_1.somaticcall",
                                  "BL23650_1_PH18_1.somaticcall",
                                  "PH29_1_T1.somaticcall",
                                  "PH29_1_PH29_T5_2_3.somaticcall",
                                  "PH29_1_PH29_T20.somaticcall",
                                  "PH29_1_T13.somaticcall",
                                  "PHD11_1_PHD11_T20.somaticcall",
                                  "PHD11_1_T7.somaticcall",
                                  "BL23650_1_T1.somaticcall",
                                  "BL23650_1_PH29_T5_2_3.somaticcall",
                                  "BL23650_1_PH29_T20.somaticcall",
                                  "BL23650_1_T13.somaticcall",
                                  "BL23650_1_PHD11_T20.somaticcall",
                                  "BL23650_1_T7.somaticcall"),
                   type= c(rep("control", 4),
                           rep("primary", 6),
                           rep("host", 6)))
comp[, cdition:= tstrsplit(comparison, "_vs_", keep= 1)]
comp[, file:= list.files("db/DNA_analysis_novogene/somatic_mutations/", 
                         paste0(file_pattern, "_", class, ".vcf.gz$"), 
                         recursive = T, 
                         full.names = T), file_pattern]
dat <- comp[,{
  fread(file,
        fill= T, 
        skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS)]
}, (comp)]

dat[cdition=="control"]


pdf("pdf/DNA/genomic_mutations.pdf", width = 14, height = 7)
par(mfrow=c(2,3))
for(class in c("InDel", "SNP"))
{
  
  files <- data.table(file= list.files("db/DNA_analysis_novogene/somatic_mutations/", 
                                       paste0("_", class, ".vcf.gz$"), 
                                       recursive = T, 
                                       full.names = T))
  files[, cdition:= gsub(paste0(".somaticcall_", class, ".vcf.gz"), "", basename(file))]
  setorderv(files, "cdition")
  
  # Primary tumors
  PH18_compare <- files[grepl("^PH18", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS)]
  }, .(file, cdition)]
  
  vl_upset_plot(split(PH18_compare$V1, PH18_compare$cdition))
  
  # PHD11 Transplantations 
  PHD11_T8_vs_primary <- files[grepl("^PHD11_1_PHD11_T8_1_2_FACS", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS)]
  }, .(file, cdition)]
  
  PHD11_T8_vs_host <- files[grepl("^BL23650_1_PHD11_T8_1_2_FACS", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS)]
  }, .(file, cdition)]
  
  PHD11_T20_vs_primary <- files[grepl("^PHD11_1_PHD11_T20", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
  }, .(file, cdition)]
  
  PHD11_T20_vs_T8 <- files[grepl("^PHD11_T20_PHD11_T8_1_2_FACS", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
  }, .(file, cdition)]
  
  PHD11_T20_vs_host <- files[grepl("^BL23650_1_PHD11_T20", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
  }, .(file, cdition)]
  
  vl_upset_plot(list(PHD11_T8_vs_primary= PHD11_T8_vs_primary$V1, 
                     PHD11_T8_vs_host= PHD11_T8_vs_host$V1,
                     PHD11_T20_vs_primary= PHD11_T20_vs_primary$V1, 
                     PHD11_T20_vs_T8= PHD11_T20_vs_host$V1,
                     PHD11_T20_vs_host= PHD11_T20_vs_host$V1))
  
  # PH29 Transplantations 
  PH29_T5_vs_primary <- files[grepl("^PH29_1_PH29_T5_2_3", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS)]
  }, .(file, cdition)]
  
  PH29_T5_vs_host <- files[grepl("^BL23650_1_PH29_T5_2_3", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS)]
  }, .(file, cdition)]
  
  PH29_T20_vs_primary <- files[grepl("^PH29_1_PH29_T20", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
  }, .(file, cdition)]
  
  PH29_T20_vs_T5 <- files[grepl("^PH29_T20_PH29_T5_2_3", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
  }, .(file, cdition)]
  
  PH29_T20_vs_host <- files[grepl("^BL23650_1_PH29_T20", cdition), {
    fread(file,
          fill= T, 
          skip = "#CHROM")[FILTER=="PASS", paste0(`#CHROM`,":",POS,":",REF,":",ALT)]
  }, .(file, cdition)]
  
  vl_upset_plot(list(PH29_T5_vs_primary= PH29_T5_vs_primary$V1, 
                     PH29_T5_vs_host= PH29_T5_vs_host$V1,
                     PH29_T20_vs_primary= PH29_T20_vs_primary$V1, 
                     PH29_T20_vs_T5= PH29_T20_vs_T5$V1,
                     PH29_T20_vs_host= PH29_T20_vs_host$V1))
}
dev.off()
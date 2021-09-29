# meta <- vl_sra_metadata(GSE = "SRA009364", SRAdb = "D:/_R_data/genomes/SRAmetadb.sqlite") # works on Mac

# m1 <- fread("Rdata/modEncode_RNASeq_metadata.txt", 
#             select = c(2,8,12,13,14,15,16,17,19,20,22,23,26,27,40,41,42,51,52,53,54,57,71))
# m2 <- as.data.table(read_excel("Rdata/modEncode_RNASeq_metadata_from_paper.xls"))
# # m2 <- m2[`Strain or Stock`=="Oregon-R-modENCODE"]
# clean <- merge(m1, m2, by.x= "sample", by.y= "GEO/SRA Sample Accession")
# 
# sel <- fread("Rdata/modEncode_RNASeq_metadata.txt", 
#             select = c(2,8,12,13,14,15,16,17,19,20,22,23,26,27,40,41,42,51,52,53,54,57,71))
# sel <- sel[study_name %in% c("modENCODE D. melanogaster Developmental Total RNA-Seq",
#                              "Developmental Time Course poly(A)  RNA Profiling in D. melanogaster",
#                              "Dissected Tissue RNAseq",
#                              "Developmental Time Course poly(A)+ RNA Profiling in D. melanogaster",
#                              "modENCODE D. melanogaster poly(A) site-seq")]
# sel <- sel[study_alias %in% c("Developmental D. melanogaster poly(A)  RNA-Seq",
#                               "Dissected Tissue RNAseq",
#                               "modENCODE D. melanogaster Mixed Embryo poly(A)  RNA-Seq",
#                               "modENCODE D. melanogaster poly(A) site-seq")]
# sel <- sel[sample_alias %in% c("D.mel. mixed embryo",
#                                "pupae, 2 days after WPP",
#                                "adult female, 5 days after eclosion",
#                                "Adult Mated Female 1 day Post-eclosion Heads",
#                                "embryos, 0-2 hr",
#                                "embryos, 6-8 hr",
#                                "embryos, 12-14 hr",
#                                "embryos, 14-16 hr",
#                                "embryos, 16-18 hr",
#                                "embryos, 22-24 hr",
#                                "adult female, 1 day after eclosion",
#                                "adult male, 1 day after eclosion",
#                                "adult male, 30 days after eclosion",
#                                "pupae, 4 days after WPP",
#                                "L3 stage larvae, clear gut PS(7-9)",
#                                "L3 CNS",
#                                "embryos,  2-4 hr",
#                                "Adult Mated Female 20 days Post-eclosion Heads",
#                                "Adult Mated Female 4 days Post-eclosion Heads",
#                                "Adult Mated Female 4 days Post-eclosion Ovaries",
#                                "Adult Mated Male 1 day Post-eclosion Heads",
#                                "Adult Mated Male 20 days Post-eclosion Heads",
#                                "Adult Mated Male 4 days Post-eclosion AccessoryGlands",
#                                "Adult Mated Male 4 days Post-eclosion Heads",
#                                "Adult Mated Male 4 days Post-eclosion Testes",
#                                "Adult Mixed Male and Female 1 day Post-eclosion Carcass",
#                                "Adult Mixed Male and Female 1 day Post-eclosion Digestive System",
#                                "Adult Mixed Male and Female 20 days Post-eclosion Carcass",
#                                "Adult Mixed Male and Female 20 days Post-eclosion Digestive System",
#                                "Adult Mixed Male and Female 4 days Post-eclosion Carcass",
#                                "Adult Mixed Male and Female 4 days Post-eclosion Digestive System",
#                                "Adult Virgin Female 1 day Post-eclosion Heads",
#                                "Adult Virgin Female 20 days Post-eclosion Heads",
#                                "Adult Virgin Female 4 days Post-eclosion Heads",
#                                "Adult Virgin Female 4 days Post-eclosion Ovaries",
#                                "L3 Carcass",
#                                "L3 Digestive System",
#                                "L3 Fat body",
#                                "L3 Imaginal Discs",
#                                "L3 Salivary Glands",
#                                "WPP+2 days CNS",
#                                "WPP+2 days Fat",
#                                "WPP Fat Body",
#                                "WPP Salivary Glands",
#                                "adult female, 30 days after eclosion",
#                                "WPP",
#                                "pupae, 12 hr after WPP",
#                                "L3 stage larvae, light blue gut PS(3-6)",
#                                "pupae, 3 days after WPP",
#                                "embryos, 20-22 hr",
#                                "L3 stage larvae, 12 hr post-molt",
#                                "adult male, 5 days after eclosion",
#                                "pupae, 24 hrs after WPP",
#                                "embryos, 4-6 hr",
#                                "L1 stage larvae",
#                                "embryos, 18-20 hr",
#                                "L2 stage larvae",
#                                "L3 stage larvae",
#                                "dark blue gut PS(1-2)",
#                                "embryos, 10-12 hr",
#                                "embryos, 8-10 hr")]
# sel[, N_rep:= .N, sample_alias]
# 
# sel[, download_file:= paste0("db/fastq/RNA_modENCODE/", run, "_", gsub(" |,|[.]|-", "_", sample_alias), ".fq.gz")]
# sel[, download_file:= gsub("__", "_", download_file)]
# sel[, {
#   if(!file.exists(download_file))
#   {
#     download.file(ftp)
#   }
# }, (sel)]


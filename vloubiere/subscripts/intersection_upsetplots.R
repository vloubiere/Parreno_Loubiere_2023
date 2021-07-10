dat <- readRDS("Rdata/DNA_novogene_clean_final.rds")

pdf("pdf/mutations_intersect_upset_plots.pdf", 
    width = 35, 
    height = 10)
par(mfrow= c(2,4))

vl_upset_plot(split(dat[CLASS=="SNP", ALLELE_UNIQ_ID], 
                    dat[CLASS=="SNP", variable]))
mtext("      SNPs - ALL")
vl_upset_plot(split(dat[CLASS=="SNP" & !(loh_from_ph18), ALLELE_UNIQ_ID], 
                    dat[CLASS=="SNP" & !(loh_from_ph18), variable]))
mtext("      SNPs - Comparison to ph18")
vl_upset_plot(split(dat[CLASS=="SNP" & !(loh_from_ph29), ALLELE_UNIQ_ID], 
                    dat[CLASS=="SNP" & !(loh_from_ph29), variable]))
mtext("      SNPs - Comparison to ph29")
vl_upset_plot(split(dat[CLASS=="SNP" & !(loh_from_phd11), ALLELE_UNIQ_ID], 
                    dat[CLASS=="SNP" & !(loh_from_phd11), variable]))
mtext("      SNPs - Comparison to phd11")


vl_upset_plot(split(dat[CLASS=="INDEL", ALLELE_UNIQ_ID], 
                    dat[CLASS=="INDEL", variable]))
mtext("      INDELs - ALL")
vl_upset_plot(split(dat[CLASS=="INDEL" & !(loh_from_ph18), ALLELE_UNIQ_ID], 
                    dat[CLASS=="INDEL" & !(loh_from_ph18), variable]))
mtext("      INDELs - Comparison to ph18")
vl_upset_plot(split(dat[CLASS=="INDEL" & !(loh_from_ph29), ALLELE_UNIQ_ID], 
                    dat[CLASS=="INDEL" & !(loh_from_ph29), variable]))
mtext("      INDELs - Comparison to ph29")
vl_upset_plot(split(dat[CLASS=="INDEL" & !(loh_from_phd11), ALLELE_UNIQ_ID], 
                    dat[CLASS=="INDEL" & !(loh_from_phd11), variable]))
mtext("      INDELs - Comparison to phd11")
dev.off()
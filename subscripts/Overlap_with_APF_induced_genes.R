setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

WPP <- as.data.table(read_excel("external_data/APF-6h_vs_APF+6h_doi10.1073pnas.1900343116_s1.xlsx"))
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
fisher.test(dat$PHD9_RECOVERY, dat$symbol %in% WPP[behavior=="WT_+6hAPF_higher", Gene.ID], alternative = "greater")
fisher.test(dat$PHD9_RECOVERY, dat$symbol %in% WPP[behavior=="WT_-6hAPF_higher", Gene.ID], alternative = "greater")

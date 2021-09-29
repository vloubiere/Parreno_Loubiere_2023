setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

dat <- data.table(file= list.files("db/DNA_analysis_novogene/somatic_mutations/", 
                                   "somaticcall_SNP.vcf", 
                                   full.names = T))
dat[, cdition:= gsub(".somaticcall_SNP.vcf.gz", "", basename(file)), file]
dat <- dat[grepl("ph29_phd11", cdition)]
dat <- dat[, fread(file, skip = "#CHROM", fill= T), (dat)]
dat <- dat[FILTER=="PASS"]
dat[, ID:= .GRP, .(`#CHROM`, POS, REF, ALT)]

.l <- copy(dat)
.l <- .l[, .(IDs= .(ID)), cdition]
res <- .l$IDs
names(res) <- .l$cdition

pdf("pdf/SNP_overlap.pdf")
vl_upset_plot(res)
dev.off()
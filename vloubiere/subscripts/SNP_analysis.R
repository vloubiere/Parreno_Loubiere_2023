dat <- data.table(file= list.files("db/vcf/", ".vcf", full.names = T))
dat[, size:= file.size(file), file]
dat <- dat[size>0]
dat[, cdition:= gsub(".vcf", "", basename(file))]
dat <- dat[, fread(file), (dat)]
setkeyv(dat, "cdition")
dat[, type:= ifelse(grepl("^INDEL", INFO), "INDEL", "SNP")]

SNP <- dat[type=="SNP"]
test <- SNP["LGL_Donor_vs_refgenome"][SNP["LGL_T5_vs_refgenome"], .N==0, .EACHI, on= c("#CHROM", "POS")]$V1
SNP["LGL_T5_vs_refgenome"][test & QUAL==40]

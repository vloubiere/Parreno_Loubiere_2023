#-------------------#
# INDELs
#-------------------#
# Import INDELs
indels <- fread("db/DNA_analysis_novogene/result/05.InDel_VarDetect/raw_variants.indel.vcf.gz", skip = "#CHROM", fill= T)
colnames(indels)[1] <- "CHROM"
indels[, ID:= .I]
# Import and merge functional annotations
idfun <- fread("db/DNA_analysis_novogene/result/05.InDel_VarDetect/raw_variants.indel.avinput.variant_function.gz", 
               fill= T, 
               select = c(1,2,3,12,13,15,16))
idexonfun <- fread("db/DNA_analysis_novogene/result/05.InDel_VarDetect/raw_variants.indel.avinput.exonic_variant_function.gz", 
                   fill= T, 
                   select = c(1,2,3,12,13,15,16))
idfun[idexonfun, c("EXON_TYPE_of_MUTATION", "EXON_ANNOTATION"):= .(`i.TYPE of MUTATION`, `i.ANNOTATION`), on= "#LINE_ID"]
# List individual alleles and merge functional annotations
all_alleles <- indels[, .(CHROM, POS, REF, ALT= unlist(strsplit(ALT, ","))), ID]
all_alleles <- cbind(all_alleles, idfun[, .(ANNO_REGION, ANNOTATION, EXON_TYPE_of_MUTATION, EXON_ANNOTATION)])
# add to indels data
indels$ANNO_REGION <- sapply(split(all_alleles$ANNO_REGION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
indels$ANNOTATION <- sapply(split(all_alleles$ANNOTATION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
indels$EXON_TYPE_of_MUTATION <- sapply(split(all_alleles$EXON_TYPE_of_MUTATION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
indels$EXON_ANNOTATION <- sapply(split(all_alleles$EXON_ANNOTATION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
indels$CLASS <- "INDEL"

#-------------------#
# SNPs
#-------------------#
# Import INDELs
snps <- fread("db/DNA_analysis_novogene/result/04.SNP_VarDetect/raw_variants.snp.vcf.gz", 
              skip = "#CHROM", 
              fill= T)
colnames(snps)[1] <- "CHROM"
snps[, ID:= .I]
# Import funcitonal annotations
snpfun <- fread("db/DNA_analysis_novogene/result/04.SNP_VarDetect/raw_variants.snp.avinput.variant_function.gz", 
                fill= T, 
                select = c(1,2,3,12,13,15,16))
snpexonfun <- fread("db/DNA_analysis_novogene/result/04.SNP_VarDetect/raw_variants.snp.avinput.exonic_variant_function.gz", 
                    fill= T, 
                    select = c(1,2,3,12,13,15,16))
# Merge function and exonic function
snpfun[snpexonfun, c("EXON_TYPE_of_MUTATION", "EXON_ANNOTATION"):= .(`i.TYPE of MUTATION`, `i.ANNOTATION`), on= "#LINE_ID"]
# List individual alleles and merge functional annotations
all_alleles <- snps[, .(CHROM, POS, REF, ALT= unlist(strsplit(ALT, ","))), ID]
all_alleles <- cbind(all_alleles, snpfun[, .(ANNO_REGION, ANNOTATION, EXON_TYPE_of_MUTATION, EXON_ANNOTATION)])
# add to indels data
snps$ANNO_REGION <- sapply(split(all_alleles$ANNO_REGION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
snps$ANNOTATION <- sapply(split(all_alleles$ANNOTATION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
snps$EXON_TYPE_of_MUTATION <- sapply(split(all_alleles$EXON_TYPE_of_MUTATION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
snps$EXON_ANNOTATION <- sapply(split(all_alleles$EXON_ANNOTATION, all_alleles$ID), function(x) paste0(na.omit(x), collapse= ","))
snps$CLASS <- "SNP"

#-------------------#
# Merge final object
#-------------------#
res <- rbind(indels, snps)
res[, ALT_UNIQ_ID:= paste0(CLASS, "_", ID)]
saveRDS(res, "Rdata/DNA_novogene_full.rds")
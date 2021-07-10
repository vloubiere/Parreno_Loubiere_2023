if(!file.exists("Rdata/DNA_novogene_full.rds"))
{
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
}

# Melt data
res <- readRDS("Rdata/DNA_novogene_full.rds")
clean <- melt(res, 
              id.vars = c("CLASS", "ID", "ALT_UNIQ_ID", "CHROM", "POS", "REF", "ALT", "QUAL", "ANNO_REGION", "ANNOTATION", "EXON_TYPE_of_MUTATION", "EXON_ANNOTATION"), 
              measure.vars = c("ph18", "ph29", "ph29_t5", "phd11", "phd11_t8"))
clean[, c("allele", "allele_counts", "total_reads"):= tstrsplit(value, ":", keep = c(1,2,3))]
clean[, total_reads:= as.numeric(ifelse(total_reads %in% c(".", "NA"), 0, total_reads))]
clean[, min_total_reads:= min(total_reads), ALT_UNIQ_ID]
# Filter low reads
clean <- clean[min_total_reads>9 & !(ALT_UNIQ_ID %in% unique(clean[allele=="./.", ALT_UNIQ_ID]))]
# Identify possible loss of heterozygosity events
#ph18
loh <- clean[variable=="ph18", .(loh_pattern= CJ(strsplit(allele, "/")[[1]], strsplit(allele, "/")[[1]])[, paste0(V1, "/", V2, collapse = "|")], .(ALT_UNIQ_ID)), allele]
loh <- loh[, .(ALT_UNIQ_ID= unlist(V2)), .(allele, loh_pattern)]
clean[loh, loh_from_ph18:= i.loh_pattern, on= "ALT_UNIQ_ID"]
clean[, loh_from_ph18:= grepl(loh_from_ph18, allele), .(allele, loh_from_ph18)]
clean[, loh_from_ph18:= as.logical(loh_from_ph18)]
#ph29
loh <- clean[variable=="ph29", .(loh_pattern= CJ(strsplit(allele, "/")[[1]], strsplit(allele, "/")[[1]])[, paste0(V1, "/", V2, collapse = "|")], .(ALT_UNIQ_ID)), allele]
loh <- loh[, .(ALT_UNIQ_ID= unlist(V2)), .(allele, loh_pattern)]
clean[loh, loh_from_ph29:= i.loh_pattern, on= "ALT_UNIQ_ID"]
clean[, loh_from_ph29:= grepl(loh_from_ph29, allele), .(allele, loh_from_ph29)]
clean[, loh_from_ph29:= as.logical(loh_from_ph29)]
#phd11
loh <- clean[variable=="phd11", .(loh_pattern= CJ(strsplit(allele, "/")[[1]], strsplit(allele, "/")[[1]])[, paste0(V1, "/", V2, collapse = "|")], .(ALT_UNIQ_ID)), allele]
loh <- loh[, .(ALT_UNIQ_ID= unlist(V2)), .(allele, loh_pattern)]
clean[loh, loh_from_phd11:= i.loh_pattern, on= "ALT_UNIQ_ID"]
clean[, loh_from_phd11:= grepl(loh_from_phd11, allele), .(allele, loh_from_phd11)]
clean[, loh_from_phd11:= as.logical(loh_from_phd11)]

# FINAL
clean <- clean[, .(CLASS, 
                   ID, 
                   ALT_UNIQ_ID, 
                   CHROM, 
                   POS, 
                   REF, 
                   ALT, 
                   QUAL, 
                   ANNO_REGION, 
                   ANNOTATION, 
                   EXON_TYPE_of_MUTATION, 
                   EXON_ANNOTATION, 
                   variable, 
                   allele, 
                   allele_counts, 
                   loh_from_ph18, 
                   loh_from_ph29, 
                   loh_from_phd11)]
clean[, ALLELE_UNIQ_ID:= paste0(ALT_UNIQ_ID, "_", allele)]
saveRDS(clean, 
        "Rdata/DNA_novogene_clean_final.rds")

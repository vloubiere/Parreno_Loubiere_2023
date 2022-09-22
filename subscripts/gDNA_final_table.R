setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

import_vcf <- function(vcf_file,
                       annotation_file,
                       exonic_function_file)
{
  # Import vcf
  vcf <- vcfR::read.vcfR(vcf_file)
  dat <- as.data.table(vcfR::getFIX(vcf))
  dat <- dat[, .(id= paste0(CHROM, ":", POS, ":", REF, "_", ALT), filter= FILTER=="PASS")]
  
  # Extract allele counts and fraction
  # vcfR::queryMETA(vcf, element = 'AD')
  # vcfR::queryMETA(vcf, element = 'AF')
  AD <- as.data.table(vcfR::extract.gt(vcf, "AD"))
  cols <- paste0(rep(names(AD), each=2), c("_ref", "_alt"))
  AD <- AD[, c(tstrsplit(TUMOR, ","), tstrsplit(NORMAL, ","))]
  AD <- AD[, lapply(.SD, as.numeric)]
  setnames(AD, cols)
  # Add phyper test
  AD[, qvalue:= phyper(TUMOR_alt, 
                       TUMOR_alt+NORMAL_alt, 
                       TUMOR_ref+NORMAL_ref, 
                       TUMOR_ref+TUMOR_alt, 
                       lower.tail = F)]
  AD[, qvalue:= p.adjust(qvalue, "fdr")]
  
  # Extract allele fraction
  AF <- as.data.table(vcfR::extract.gt(vcf, "AF"))
  setnames(AF, paste0(names(AF), "_alt_freq"))
  AF <- AF[, lapply(.SD, as.numeric)]
  
  # Extract annotation
  annot <- fread(annotation_file, col.names = c("annotation", "FBtr"), sel= 1:2, fill = T)
  
  # Extract funciton
  type <- fread(exonic_function_file, sel= 1:2, col.names = c("line", "type"), fill= T)
  
  # RETURN
  dat <- cbind(dat, AD, AF, annot)
  dat[as.numeric(gsub("line", "", type$line)), type:= type$type]
  return(dat[(filter), !"filter"])
}
dat <- data.table(cdition= list.files("db/DNA_analysis_novogene/", 
                                      ".somaticcall_SNP.vcf.gz$", 
                                      recursive = T))
dat[, pattern:= tstrsplit(cdition, ".somaticcall", keep= 1), cdition]
dat[, cdition:= gsub("PH18_1_", "", pattern)]
dat <- dat[, .(alt_class= c("SNP", "InDel")), (dat)]
dat <- dat[, import_vcf(list.files("db/DNA_analysis_novogene/", paste0(pattern, ".*", alt_class, ".vcf.gz"), full.names = T),
                        list.files("db/DNA_analysis_novogene/", paste0(pattern, ".*", alt_class, ".avinput.variant_function.gz"), full.names = T),
                        list.files("db/DNA_analysis_novogene/", paste0(pattern, ".*", alt_class, ".avinput.exonic_variant_function.gz"), full.names = T)), (dat)]
# PH18
dat[, PH18:= any(grepl("PH18", cdition)), id]
# Define tumor specific vs enriched
dat[, class:= cut(NORMAL_alt_freq, 
                  c(-Inf, 0, Inf), 
                  c("Tumor specific", "Tumor enriched"))]
dat[, col:= c("darkorange2", "dodgerblue1")[class]]
dat[!(PH18), occurence:= ifelse(.N>1, "shared among several conditions", "single condition"), id]

# SAVE
fwrite(dat,
       "Rdata/gDNA_final_table.txt",
       sep= "\t",
       na= NA)

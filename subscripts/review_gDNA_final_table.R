setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(data.table)

import_vcf <- function(vcf_file,
                       annotation_file,
                       exonic_function_file)
{
  # Import vcf
  vcf <- vcfR::read.vcfR(vcf_file)
  dat <- as.data.table(vcfR::getFIX(vcf))
  dat <- dat[, .(id= paste0("chr", CHROM, ":", POS, "-", as.numeric(POS)+nchar(REF), ":", REF, "_", ALT), 
                 filter= FILTER=="PASS")]

  # Extract allele counts ----
  AD <- as.data.table(vcfR::extract.gt(vcf, "AD"))
  cols <- paste0(rep(names(AD), each=2), c("_ref", "_alt"))
  AD <- AD[, c(tstrsplit(TUMOR, ","), tstrsplit(NORMAL, ","))]
  AD <- AD[, lapply(.SD, as.numeric)]
  setnames(AD, cols)
  dat <- cbind(dat, AD)
  
  # Extract annotation ----
  annot <- fread(annotation_file,
                 sel= c(1,2,11,12,14,15),
                 fill = T)
  annot[, id:= paste0("chr", CHROM, ":", POS, "-", as.numeric(POS)+nchar(REF), ":", REF, "_", ALT)]
  annot[, ANNOTATION:= lapply(strsplit(ANNOTATION, ",|;"), function(x) sapply(strsplit(x, "\\("), "[", 1))]
  dat <- merge(dat,
               annot[, .(id, 
                         annotation= `#ANNO_REGION`,
                         FBtr= ANNOTATION)],
               by= "id",
               all.x= T)
  
  # Extract function ----
  type <- fread(exonic_function_file,
                sel= c(2,12,13,15,16),
                fill= T)
  type[, id:= paste0("chr", CHROM, ":", POS, "-", as.numeric(POS)+nchar(REF), ":", REF, "_", ALT)]
  dat <- merge(dat,
               type[, .(id,
                        exonic_annotation= `TYPE of MUTATION`)],
               by= "id",
               all.x= T)
  
  # RETURN ----
  return(dat)
}

# Import data ----
meta <- data.table(vcf_file= list.files("db/DNA_analysis_novogene/review/",
                                        "somaticcall_SNP.vcf.gz|somaticcall_InDel.vcf.gz", recursive = T, full.names = T))
meta[, annotation_file:= gsub(".vcf.gz$", ".avinput.variant_function.gz", vcf_file)]
meta[, exonic_function_file:= gsub(".vcf.gz$", ".avinput.exonic_variant_function.gz", vcf_file)]
meta[, cdition:= tstrsplit(basename(vcf_file), ".somaticcall", keep= 1), vcf_file]
meta[, ctl:= fcase(grepl("^PH18_1", cdition), "PH18_1",
                  grepl("^PH18_5", cdition), "PH18_5")]
meta[, cdition:= gsub(paste0(ctl, "_"), "", cdition), ctl]
meta[, alt_class:= fcase(grepl("_SNP", vcf_file), "SNP",
                        grepl("_InDel", vcf_file), "InDel")]

# Check if the two samples are from the same batch ----
meta[, matching_batch:= 
       (grepl("^PH18_2$|^PH29|^PHD11|^PHD9", cdition) & ctl=="PH18_1")
     | (grepl("^PH18_6$|^L3_E_", cdition) & ctl=="PH18_5")]
meta <- meta[(matching_batch), !"matching_batch"]

# Import vcf files ----
dat <- meta[, {
  import_vcf(vcf_file, annotation_file, exonic_function_file)
}, .(cdition, ctl, alt_class)]

# Remove mutations that don't pass the filter or from non-canonical chr ----
dat <- dat[(filter), !"filter"]
dat <- dat[grepl("^chr2L:|^chr2R:|^chr3L:|^chr3R:|^chr4:|^chrX:|^chrY:", id)]

# Retrieve gene symbols ----
transcripts <- rtracklayer::import("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(transcripts) <- "UCSC"
transcripts <- as.data.table(transcripts)
transcripts <- transcripts[!is.na(transcript_id), .(FBtr= transcript_id, symbol= gene_symbol)]
transcripts <- unique(transcripts)
exonic_mut <- dat[annotation=="exonic" & exonic_annotation!="synonymous SNV", .(FBtr= unlist(FBtr)), id]
exonic_mut <- unique(exonic_mut)
exonic_mut <- merge(exonic_mut,
                    transcripts)
setorderv(exonic_mut, "symbol")
exonic_mut <- exonic_mut[, .(symbol= paste0(unique(symbol), collapse= ",")), id]
dat <- merge(dat,
             exonic_mut,
             all.x= T)

# Compute allele frequency ----
dat[, allele_freq:= TUMOR_alt/(TUMOR_alt+TUMOR_ref)]

# Compute number of tumors with high allelic freq
dat[, nCtl:= sum(grepl("^PH18", cdition) & allele_freq>0.2), id]
dat[, nTum:= sum(allele_freq>0.2)-nCtl, id]
tot <- length(unique(grep("^PH18", dat$cdition, value = T, invert = T)))
dat[nCtl==0 & nTum>0, tumor_fraction:= factor(paste0(nTum, "/", tot), paste0(seq(12), "/", tot))]

# SAVE ----
dat <- dat[, .(id, cdition, allele_freq, tumor_fraction, annotation, symbol)]
dat[, cdition:= gsub("^(.*)_(.*)$", "\\1.\\2", unique(cdition)), cdition]
dat[, cdition:= gsub("PH18", "no ph-KD", unique(cdition)), cdition]
dat[, cdition:= gsub("PH29", "Constant ph-KD", unique(cdition)), cdition]
dat[, cdition:= gsub("PHD9", "Trans. ph-KD d9", unique(cdition)), cdition]
dat[, cdition:= gsub("PHD11", "Trans. ph-KD d11", unique(cdition)), cdition]
dat[, cdition:= gsub("L3_E_24", "eL3 1d recovery", unique(cdition)), cdition]
dat[, cdition:= gsub("L3_E_96", "eL3 4d recovery", unique(cdition)), cdition]
dat[, cdition:= gsub("L3_E_144", "eL3 6d recovery", unique(cdition)), cdition]
dat[, cdition:= factor(cdition,
                       c("no ph-KD.2",
                         "Constant ph-KD.1",
                         "Constant ph-KD.2",
                         "Trans. ph-KD d9.1",
                         "Trans. ph-KD d9.2",
                         "Trans. ph-KD d11.1",
                         "Trans. ph-KD d11.2",
                         "no ph-KD.6",
                         "eL3 1d recovery.1",
                         "eL3 1d recovery.2",
                         "eL3 4d recovery.1",
                         "eL3 4d recovery.2",
                         "eL3 6d recovery.1",
                         "eL3 6d recovery.2"))]
saveRDS(dat,
        "Rdata/review_gDNA_final_table.rds")

setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- data.table(vcf_file= list.files("db/DNA_analysis_novogene/review/07.CNV_VarDetect/", ".vcf.gz", recursive = T, full.names = T))
dat[, annot_file:= list.files("db/DNA_analysis_novogene/review/07.CNV_VarDetect/", ".cnv.avinput.variant_function.gz", recursive = T, full.names = T)]
dat[, cdition:= tstrsplit(basename(vcf_file), ".cnv", keep= 1), vcf_file]
dat <- dat[, {
  vcf <- vcfR::read.vcfR(vcf_file)
  .c1 <- as.data.table(vcfR::getFIX(vcf))
  .c1[, `#LINE_ID`:= paste0("line", .I)]
  .c2 <- fread(annot_file,
               fill= T)
  .c <- merge(.c1,
              .c2,
              by= "#LINE_ID")
  .c[, id:= paste0(CHROM.x, ":", POS1, "-", POS2, "_", `TYPE of CNV`)]
  .c[, ANNOTATION:= lapply(strsplit(ANNOTATION, ",|;"), function(x) sapply(strsplit(x, "\\("), "[", 1))]
  .c[, .(id, `TYPE of CNV`, NORMALIZED_RD, ANNO_REGION, ANNOTATION)]
}, cdition]

# Select canonical chromosomes ----
dat <- dat[grepl("^2L:|^2R:|^3L:|^3R:|^4:|^X:|^Y:", id)]

# Annotation ----
dat[, annotation:= `TYPE of CNV`]

# Normalized ratio
dat[, norm:= tstrsplit(NORMALIZED_RD, ":", keep= 2, type.convert = T)]

# Retrieve gene symbols ----
transcripts <- rtracklayer::import("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(transcripts) <- "UCSC"
transcripts <- as.data.table(transcripts)
transcripts <- transcripts[!is.na(transcript_id), .(FBtr= transcript_id, symbol= gene_symbol)]
transcripts <- unique(transcripts)
genes <- dat[lengths(ANNOTATION)>0, .(FBtr= unlist(ANNOTATION)), id]
genes <- unique(genes)
genes <- merge(genes,
               transcripts)
setorderv(genes, "symbol")
genes <- genes[, .(symbol= paste0(unique(symbol), collapse= ",")), id]
dat <- merge(dat,
             genes,
             all.x= T)

# Compute fraction of tumors ----
dat[, nCtl:= sum(grepl("^PH18", cdition)), id]
dat[, nTum:= .N-nCtl, id]
tot <- length(unique(grep("^PH18", dat$cdition, value = T, invert = T)))
dat[nCtl==0, tumor_fraction:= factor(paste0(nTum, "/", tot), paste0(seq(12), "/", tot))]

# SAVE ----
dat <- dat[, .(id, cdition, norm, tumor_fraction, annotation, ANNO_REGION, symbol)]
dat[, cdition:= gsub("^(.*)_(.*)$", "\\1.\\2", unique(cdition)), cdition]
dat[, cdition:= gsub("PH18", "no ph-KD", unique(cdition)), cdition]
dat[, cdition:= gsub("PH29", "Constant ph-KD", unique(cdition)), cdition]
dat[, cdition:= gsub("PHD9", "Trans. ph-KD d9", unique(cdition)), cdition]
dat[, cdition:= gsub("PHD11", "Trans. ph-KD d11", unique(cdition)), cdition]
dat[, cdition:= gsub("L3_E_24", "eL3 1d recovery", unique(cdition)), cdition]
dat[, cdition:= gsub("L3_E_96", "eL3 4d recovery", unique(cdition)), cdition]
dat[, cdition:= gsub("L3_E_144", "eL3 6d recovery", unique(cdition)), cdition]
dat[, cdition:= factor(cdition,
                       c("no ph-KD.1",
                         "no ph-KD.2",
                         "no ph-KD.5",
                         "no ph-KD.6",
                         "Constant ph-KD.1",
                         "Constant ph-KD.2",
                         "Trans. ph-KD d9.1",
                         "Trans. ph-KD d9.2",
                         "Trans. ph-KD d11.1",
                         "Trans. ph-KD d11.2",
                         "eL3 1d recovery.1",
                         "eL3 1d recovery.2",
                         "eL3 4d recovery.1",
                         "eL3 4d recovery.2",
                         "eL3 6d recovery.1",
                         "eL3 6d recovery.2"))]
saveRDS(dat,
        "Rdata/review_gDNA_CNVs_final_table.rds")

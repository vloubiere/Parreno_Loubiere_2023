
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

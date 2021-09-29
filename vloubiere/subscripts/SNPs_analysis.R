#----------------------------#
# SNPs
#----------------------------#
snp <- fread("db/DNA_analysis_novogene/result/04.SNP_VarDetect/raw_variants.snp.vcf.gz", skip = "#CHROM", fill= T)
snp[, line_number:= .I]
# Differences
snpd <- melt(snp, id.vars = "line_number", measure.vars = patterns("ph"))
snpd[, alleles:= tstrsplit(value, ":", keep = 1)]
snpd[, total_reads:= tstrsplit(value, ":", keep = 3)]
snpd[total_reads==".", total_reads:= NA]
snpd[, total_reads:= as.numeric(total_reads)]
snpd[, min_total_reads:= min(total_reads), line_number]
snpd[, ID:= paste0(line_number, "_", alleles)]
# Filter low reads
snpd <- snpd[min_total_reads>9]
# Loss of heterozygosity
loh <- snpd[variable=="ph18", .(loh_pattern= CJ(strsplit(alleles, "/")[[1]], strsplit(alleles, "/")[[1]])[, paste0(V1, "/", V2, collapse = "|")], .(line_number)), alleles]
loh <- loh[, .(line_number= unlist(V2)), .(alleles, loh_pattern)]
snpd[loh, loh_pattern:= i.loh_pattern, on= "line_number"]
snpd[, loh_check:= !grepl(alleles, loh_pattern), .(alleles, loh_pattern)]

snpall <- lapply(split(snpd, snpd$variable), function(x) x$ID)
snpspecific <- lapply(split(snpd[(loh_check)], snpd[(loh_check),variable]), function(x) x$ID)

# Overlaps
pdf("pdf/snp_intersection.pdf", width = 10, height = 10)
par(mfrow=c(2,1))
vl_upset_plot(snpall)
vl_upset_plot(snpspecific[-1])
dev.off()

# snp variant function
snpfun <- fread("db/DNA_analysis_novogene/result/04.SNP_VarDetect/raw_variants.snp.avinput.variant_function.gz", fill= T)
grid <- data.table(type= unique(snpfun[, ANNO_REGION]), 
                   Cc= matlab.like2(length(unique(snpfun[, ANNO_REGION]))))

pdf("pdf/pie_charts_SNPs_function.pdf", height = 10)
par(mfrow=c(3,2))
.c <- snpfun[grid, .(.N, Cc), .EACHI, on= "ANNO_REGION==type"]
pie(.c$N, 
    labels = .c$ANNO_REGION,
    col = .c$Cc,
    main= paste0("All SNPs (", formatC(sum(.c$N), big.mark = ","), ")"))
ph18_snp <- paste0("line", snpd[variable=="ph18" & !(alleles %in% c("0/0", "./.")), line_number])
.c <- snpfun[`#LINE_ID` %in% ph18_snp][grid, .(.N, Cc), .EACHI, on= "ANNO_REGION==type"]
pie(.c$N, 
    labels = .c$ANNO_REGION,
    col = .c$Cc,
    main= paste0("ph18 SNPs (", formatC(sum(.c$N), big.mark = ","), ")"))
for(cdition in c("ph29", "phd11", "ph29_t5", "phd11_t8"))
{
  sel <- paste0("line", gsub("(.*)_.*", "\\1", snpspecific[[cdition]]))
  .c <- snpfun[`#LINE_ID` %in% sel][grid, .(.N, Cc), .EACHI, on= "ANNO_REGION==type"]
  pie(.c$N, 
      labels = .c$ANNO_REGION,
      col = .c$Cc,
      main= paste0(cdition, "-specific SNPs (", formatC(sum(.c$N), big.mark = ","), ")"))
}
dev.off()

# snp exon variant function
snpexonfun <- fread("db/DNA_analysis_novogene/result/04.SNP_VarDetect/raw_variants.snp.avinput.exonic_variant_function.gz", fill= T)
grid <- data.table(type= unique(snpexonfun[, `TYPE of MUTATION`]), 
                   Cc= matlab.like2(length(unique(snpexonfun[, `TYPE of MUTATION`]))))

pdf("pdf/pie_charts_SNPs_exon_function.pdf", height = 10)
par(mfrow=c(3,2))
.c <- snpexonfun[grid, .(.N, Cc), .EACHI, on= "TYPE of MUTATION==type"]
pie(.c$N, 
    labels = .c$`TYPE of MUTATION`,
    col = .c$Cc,
    main= paste0("All SNPs (", formatC(sum(.c$N), big.mark = ","), ")"))
ph18_snp <- paste0("line", snpd[variable=="ph18" & !(alleles %in% c("0/0", "./.")), line_number])
.c <- snpexonfun[`#LINE_ID` %in% ph18_snp][grid, .(.N, Cc), .EACHI, on= "TYPE of MUTATION==type"]
pie(.c$N, 
    labels = .c$`TYPE of MUTATION`,
    col = .c$Cc,
    main= paste0("ph18 SNPs (", formatC(sum(.c$N), big.mark = ","), ")"))
for(cdition in c("ph29", "phd11", "ph29_t5", "phd11_t8"))
{
  sel <- paste0("line", gsub("(.*)_.*", "\\1", snpspecific[[cdition]]))
  .c <- snpexonfun[`#LINE_ID` %in% sel][grid, .(.N, Cc), .EACHI, on= "TYPE of MUTATION==type"]
  pie(.c$N, 
      labels = .c$`TYPE of MUTATION`,
      col = .c$Cc,
      main= paste0(cdition, "-specific SNPs (", formatC(sum(.c$N), big.mark = ","), ")"))
}
dev.off()

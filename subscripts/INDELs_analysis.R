#----------------------------#
# INDELS
#----------------------------#
indels <- fread("db/DNA_analysis_novogene/result/05.InDel_VarDetect/raw_variants.indel.vcf.gz", skip = "#CHROM", fill= T)
indels[, line_number:= .I]

# Differences
iddiff <- melt(indels, id.vars = "line_number", measure.vars = patterns("ph"))
iddiff[, alleles:= tstrsplit(value, ":", keep = 1)]
iddiff[, total_reads:= tstrsplit(value, ":", keep = 3)]
iddiff[total_reads==".", total_reads:= NA]
iddiff[, total_reads:= as.numeric(total_reads)]
iddiff[, min_total_reads:= min(total_reads), line_number]
iddiff[, ID:= paste0(line_number, "_", alleles)]
# Filter low reads
iddiff <- iddiff[min_total_reads>9]
# Loss of heterozygosity
loh <- iddiff[variable=="ph18", .(loh_pattern= CJ(strsplit(alleles, "/")[[1]], strsplit(alleles, "/")[[1]])[, paste0(V1, "/", V2, collapse = "|")], .(line_number)), alleles]
loh <- loh[, .(line_number= unlist(V2)), .(alleles, loh_pattern)]
iddiff[loh, loh_pattern:= i.loh_pattern, on= "line_number"]
iddiff[, loh_check:= !grepl(alleles, loh_pattern), .(alleles, loh_pattern)]

idall <- lapply(split(iddiff, iddiff$variable), function(x) x$ID)
idspecific <- lapply(split(iddiff[(loh_check)], iddiff[(loh_check),variable]), function(x) x$ID)

pdf("pdf/indel_intersection.pdf", width = 10, height = 10)
par(mfrow=c(2,1))
vl_upset_plot(idall)
vl_upset_plot(idspecific[-1])
dev.off()

# Indel variant function
idfun <- fread("db/DNA_analysis_novogene/result/05.InDel_VarDetect/raw_variants.indel.avinput.variant_function.gz", fill= T)
grid <- data.table(type= unique(idfun[, ANNO_REGION]), 
                   Cc= matlab.like2(length(unique(idfun[, ANNO_REGION]))))

pdf("pdf/pie_charts_INDELs_function.pdf", height = 10)
par(mfrow=c(3,2))
.c <- idfun[grid, .(.N, Cc), .EACHI, on= "ANNO_REGION==type"]
pie(.c$N, 
    labels = .c$ANNO_REGION,
    col = .c$Cc,
    main= paste0("All INDELs (", formatC(sum(.c$N), big.mark = ","), ")"))
ph18_sel <-  paste0("line", iddiff[variable=="ph18" & !(alleles %in% c("0/0", "./.")), line_number])
.c <- idfun[`#LINE_ID` %in% ph18_sel][grid, .(.N, Cc), .EACHI, on= "ANNO_REGION==type"]
pie(.c$N, 
    labels = .c$ANNO_REGION,
    col = .c$Cc,
    main= paste0("ph18 INDELs (", formatC(sum(.c$N), big.mark = ","), ")"))
for(cdition in c("ph29", "phd11", "ph29_t5", "phd11_t8"))
{
  sel <- paste0("line", gsub("(.*)_.*", "\\1", idspecific[[cdition]]))
  .c <- idfun[`#LINE_ID` %in% sel][grid, .(.N, Cc), .EACHI, on= "ANNO_REGION==type"]
  pie(.c$N, 
      labels = .c$ANNO_REGION,
      col = .c$Cc,
      main= paste0(cdition, "-specific INDELs (", formatC(sum(.c$N), big.mark = ","), ")"))
}
dev.off()

# Indel exon variant function
idexonfun <- fread("db/DNA_analysis_novogene/result/05.InDel_VarDetect/raw_variants.indel.avinput.exonic_variant_function.gz", fill= T)
grid <- data.table(type= unique(idexonfun[, `TYPE of MUTATION`]), 
                   Cc= matlab.like2(length(unique(idexonfun[, `TYPE of MUTATION`]))))

pdf("pdf/pie_charts_INDELs_exon_function.pdf", height = 10)
par(mfrow=c(3,2))
.c <- idexonfun[grid, .(.N, Cc), .EACHI, on= "TYPE of MUTATION==type"]
pie(.c$N, 
    labels = .c$`TYPE of MUTATION`,
    col = .c$Cc,
    main= paste0("All INDELs (", formatC(sum(.c$N), big.mark = ","), ")"))
ph18_sel <-  paste0("line", iddiff[variable=="ph18" & !(alleles %in% c("0/0", "./.")), line_number])
.c <- idexonfun[`#LINE_ID` %in% ph18_sel][grid, .(.N, Cc), .EACHI, on= "TYPE of MUTATION==type"]
pie(.c$N, 
    labels = .c$`TYPE of MUTATION`,
    col = .c$Cc,
    main= paste0("ph18 INDELs (", formatC(sum(.c$N), big.mark = ","), ")"))
for(cdition in c("ph29", "phd11", "ph29_t5", "phd11_t8"))
{
  sel <- paste0("line", gsub("(.*)_.*", "\\1", idspecific[[cdition]]))
  .c <- idexonfun[`#LINE_ID` %in% sel][grid, .(.N, Cc), .EACHI, on= "TYPE of MUTATION==type"]
  pie(.c$N, 
      labels = .c$`TYPE of MUTATION`,
      col = .c$Cc,
      main= paste0(cdition, "-specific INDELs (", formatC(sum(.c$N), big.mark = ","), ")"))
}
dev.off()




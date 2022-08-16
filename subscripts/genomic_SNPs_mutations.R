setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)
require(readxl)

if(!exists("dat"))
{
  # Collect files paths
  dat <- data.table(pattern= c("PH18.*PH18.*.somaticcall_SNP.vcf.gz$",
                               "PH18.*PH29.*.somaticcall_SNP.vcf.gz$",
                               "PH18.*PHD9.*.somaticcall_SNP.vcf.gz$",
                               "PH18.*PHD11.*.somaticcall_SNP.vcf.gz$"))
  dat <- dat[, .(file= list.files("db/DNA_analysis_novogene/somatic_mutations/", 
                                  pattern,
                                  recursive = T,
                                  full.names = T)), pattern]
  dat[, annot:= gsub(".vcf.gz$",
                     ".avinput.variant_function.gz",
                     file), file]
  dat[, exon_annot:= gsub(".vcf.gz$",
                          ".avinput.exonic_variant_function.gz",
                          file), file]
  dat[, cdition:= gsub(".somaticcall_SNP.vcf.gz$", "", basename(file))]
  # Import files
  dat <- dat[, {
    # Import original vcf file
    .c <- fread(file,
                fill= T, 
                skip = "#CHROM", 
                sel= c(1,2,4,5,7,9,10,11),
                col.names = c("CHROM", "POS", "REF", "ALT", "FILTER", "FORMAT", "TUMOR", "NORMAL"))
    # Only keep somatic mutations passing all filters
    .c <- .c[FILTER=="PASS"]
    # Reshape NORM and TUM columns before binding
    NORM <- mcmapply(function(names, var) {
      x <- as.data.table(tstrsplit(var, ":"))
      setnames(x, paste0(unlist(tstrsplit(names, ":")), "_NORM"))
      return(x)
    }, 
    names= .c$FORMAT, 
    var= .c$NORMAL, 
    SIMPLIFY = F)
    TUM <- mcmapply(function(names, var) {
      x <- as.data.table(tstrsplit(var, ":"))
      setnames(x, paste0(unlist(tstrsplit(names, ":")), "_TUM"))
      return(x)
    }, 
    names= .c$FORMAT, 
    var= .c$TUMOR, 
    SIMPLIFY = F)
    # Make object
    .c <- cbind(.c[, CHROM:ALT], 
                rbindlist(NORM, fill= T), 
                rbindlist(TUM, fill= T))
    # Add functional annotations
    .a <- fread(annot, 
                fill= T,
                sel= c(11,12,14,15,1,2),
                col.names = c("CHROM", "POS", "REF", "ALT", "ANNOT_REGION", "ANNOTATION"))
    .c <- merge(.c, 
                .a,
                by= c("CHROM", "POS", "ALT"), 
                all.x= T)
    # Add Exonic annotation
    .e <- fread(exon_annot, 
                fill= T,
                sel= c(12,13,15,16,2,3),
                col.names = c("CHROM", "POS", "REF", "ALT", "MUTATION_TYPE", "EXON_ANNOTATION"))
    .c <- merge(.c, 
                .e,
                by= c("CHROM", "POS", "ALT"),
                all.x= T)
  }, (dat)]
  dat[, ID:= paste0(CHROM, "_", POS, "_", ALT)]
  # Remove all mutations which are found by comparing the two 18C conditions
  dat <- dat[!ID %in% dat[cdition=="PH18_1_PH18_2", ID]]
}

# SNPs overlaps
pdf("pdf/DNA/gDNA_SNPs_overlaps.pdf", 
    width = 20, 
    height = 7)
vl_upset_plot(split(dat$ID, dat$cdition))
title("all SNPs")
vl_upset_plot(split(dat[MUTATION_TYPE %in% c("nonsynonymous SNV", "stopgain"), ID], 
                    dat[MUTATION_TYPE %in% c("nonsynonymous SNV", "stopgain"), cdition]))
title("non sysnonymous/stopgain mutations only")
dev.off()

# SNPs overlaps
cols <- c("AF_NORM", "AF_TUM")
dat[, (cols):= lapply(.SD, as.numeric), .SDcols= cols]
pdf("pdf/DNA/gDNA_SNPs_Allele_freq.pdf", 
    width = 14, 
    height = 7)
par(mfrow=c(2,3))
dat[, {
  plot(AF_NORM, 
       AF_TUM, 
       main= cdition,
       ) 
}, cdition]
dev.off()

# SNPs functional features
counts <- dcast(dat, ANNOT_REGION~cdition, fun.aggregate = length)
counts <- as.matrix(counts, 1)
counts <- counts[order(apply(counts, 1, sum), decreasing = T),]
counts.e <- dcast(dat[!is.na(MUTATION_TYPE)], 
                  MUTATION_TYPE~cdition, 
                  fun.aggregate = length)
counts.e <- as.matrix(counts.e, 1)
counts.e <- counts.e[order(apply(counts.e, 1, sum), decreasing = T),]


pdf("pdf/DNA/gDNA_SNPs_annotation.pdf")
par(las= 2,
    mar= c(9,4,2,2))
barplot(counts, 
        col= rev(vl_palette_many_categ(nrow(counts))),
        ylab= "SNP counts")
legend("topleft",
       bty= "n",
       rownames(counts),
       fill= rev(vl_palette_many_categ(nrow(counts))))
barplot(counts.e, 
        col= rev(vl_palette_many_categ(nrow(counts.e))),
        ylab= "Exonic SNP counts")
legend("topleft",
       bty= "n",
       rownames(counts.e),
       fill= rev(vl_palette_many_categ(nrow(counts.e))))
dev.off()

# Counts mutations per gene
gene_counts <- dat[MUTATION_TYPE %in% c("nonsynonymous SNV", "stopgain"), .(cdition, CHROM, POS, ALT, ANNOTATION, MUTATION_TYPE)]
gene_counts <- gene_counts[, .(gene= strsplit(ANNOTATION, ",")[[1]]), gene_counts]
FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
FBGN <- as.data.table(FBGN)
FBGN <- FBGN[, .(gene_id, transcript_id, gene_symbol)]
gene_counts[FBGN, c("FBgn", "symbol"):= .(i.gene_id, i.gene_symbol), on= "gene==transcript_id"]
gene_counts <- unique(gene_counts[!is.na(FBgn), .(cdition, CHROM, POS, ALT, MUTATION_TYPE, FBgn, symbol)])
gene_counts <- dcast(gene_counts, 
                     FBgn+symbol~cdition,
                     fun.aggregate = length)
cols <- grep("^PH18", names(gene_counts), value = T)
gene_counts[, total:= rowSums(.SD), .SDcols= cols]
setorderv(gene_counts, "total", -1)
fwrite(gene_counts, 
       "db/DNA_analysis_novogene/table_SNPs_nonsynonymous_stopgain_mutations.txt")

pdf("pdf/DNA/gDNA_function_annotation_SNPs.pdf")
par(mar= c(5,20,5,10))
vl_GO_enrich(gene_counts[total>0, FBgn])
# vl_GO_enrich(gene_counts[total>1, FBgn])
# vl_GO_enrich(gene_counts[PH18_1_PH29_1>0 | PH18_1_PH29_2>0, FBgn])
# test <- vl_GO_enrich(gene_counts[PH18_1_PHD11_1>0 | PH18_1_PHD11_2>0, FBgn], plot=F)
# vl_GO_enrich(gene_counts[PH18_1_PHD9_1>0 | PH18_1_PHD9_2>0, FBgn])
dev.off()

# Overlap between mutated genes
genes_mut <- unlist(as.list(gene_counts[, lapply(.SD, function(x) list(FBgn[x>0])), .SDcols= patterns("^PH18")]), recursive = F)

pdf("pdf/DNA/gDNA_mutated_genes_overlaps_SNPs.pdf",
    width= 18)
vl_upset_plot(genes_mut)
dev.off()



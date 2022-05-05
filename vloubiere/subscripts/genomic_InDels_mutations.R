setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)
require(readxl)

# Import data
dat <- data.table(pattern= c("PH18.*PH18.*.somaticcall_InDel.vcf.gz$",
                             "PH18.*PH29.*.somaticcall_InDel.vcf.gz$",
                             "PH18.*PHD9.*.somaticcall_InDel.vcf.gz$",
                             "PH18.*PHD11.*.somaticcall_InDel.vcf.gz$"))
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
dat[, cdition:= gsub(".somaticcall_InDel.vcf.gz$", "", basename(file))]
dat <- dat[, {
  .c <- fread(file,
              fill= T, 
              skip = "#CHROM", 
              sel= c(1,2,5,7),
              col.names = c("CHROM", "POS", "ALT", "FILTER"))
  .c <- .c[FILTER=="PASS"]
  .a <- fread(annot, 
              fill= T,
              sel= c(11,12,15,1,2),
              col.names = c("CHROM", "POS", "ALT", "ANNOT_REGION", "ANNOTATION"))
  .c <- merge(.c, 
              .a,
              by= c("CHROM", "POS", "ALT"), 
              all.x= T)
  .e <- fread(exon_annot, 
              fill= T,
              sel= c(12,13,16,2,3),
              col.names = c("CHROM", "POS", "ALT", "MUTATION_TYPE", "EXON_ANNOTATION"))
  .c <- merge(.c, 
              .e,
              by= c("CHROM", "POS", "ALT"),
              all.x= T)
}, (dat)]
dat[, ID:= paste0(CHROM, "_", POS, "_", ALT)]
dat <- dat[!ID %in% dat[cdition=="PH18_1_PH18_2", ID]]

# SNPs overlaps
pdf("pdf/DNA/gDNA_InDels_overlaps.pdf", 
    width = 20, 
    height = 7)
vl_upset_plot(split(dat$ID, dat$cdition))
title("all InDels")
vl_upset_plot(split(dat[!is.na(MUTATION_TYPE), ID], dat[!is.na(MUTATION_TYPE), cdition]))
title("Pontentially deleterious InDels")
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


pdf("pdf/DNA/gDNA_InDels_annotation.pdf")
par(las= 2,
    mar= c(9,4,2,2))
barplot(counts, 
        col= rev(vl_palette_many_categ(nrow(counts))),
        ylab= "InDel counts")
legend("topleft",
       bty= "n",
       rownames(counts),
       fill= rev(vl_palette_many_categ(nrow(counts))))
barplot(counts.e, 
        col= rev(vl_palette_many_categ(nrow(counts.e))),
        ylab= "Exonic InDel counts")
legend("topleft",
       bty= "n",
       rownames(counts.e),
       fill= rev(vl_palette_many_categ(nrow(counts.e))))
dev.off()

# Counts mutations per gene
gene_counts <- dat[!is.na(MUTATION_TYPE) , .(cdition, CHROM, POS, ALT, ANNOTATION, MUTATION_TYPE)]
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
       "db/DNA_analysis_novogene/table_InDels_pontentially_deleterious_mutations.txt")

GOs <- vl_GO_enrich(gene_counts[total>0, FBgn], 
                    plot= F)
pdf("pdf/DNA/gDNA_function_annotation_InDels.pdf")
par(mar= c(5,20,2,10))
plot(GOs)
dev.off()

# Overlap between mutated genes
genes_mut <- unlist(as.list(gene_counts[, lapply(.SD, function(x) list(FBgn[x>0])), .SDcols= patterns("^PH18")]), recursive = F)

pdf("pdf/DNA/gDNA_mutated_genes_overlaps_InDels.pdf",
    width= 18)
vl_upset_plot(genes_mut)
dev.off()

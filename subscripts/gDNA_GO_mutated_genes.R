# setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)
require(vlfunctions)

######################################################################
# Import data
######################################################################
genes <- rtracklayer::import("../../genomes/flybase/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)

dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[(filter) & exonic_annotation %in% c("frameshift deletion", "frameshift insertion", "nonsynonymous SNV", "stopgain")
           & cdition !="PH18_2"]

FBgn <- unique(unlist(dat[, FBgn]))
enr1 <- vl_GO_enrich(list(FBgn= FBgn), species = "Dm")

FBgn <- dat[, .(FBgn= unique(unlist(FBgn))), cdition]
FBgn <- unique(FBgn[, .(.N>=2), FBgn][(V1), FBgn])
enr2 <- vl_GO_enrich(list(FBgn= FBgn), species = "Dm")

pdf("pdf/gDNA_mutant_genes_GO.pdf", 4, 4)
par(mar= c(7,25,3,7),
    las= 2,
    cex= 0.5)
pl <- plot(enr1, padj_cutoff= 0.05)
pl <- plot(enr2, padj_cutoff= 0.05)
dev.off()

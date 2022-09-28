setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
genes <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)

dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!(PH18) & type %in% c("nonsynonymous SNV", "stopgain")]
dat <- dat[, .(FBgn= unlist(FBgn)), .(id, class, occurence, cdition)]
dat <- na.omit(dat)
dat[, occurence:= factor(occurence,
                         c("single condition",
                           "shared >=1 conditions"))]

pdf("pdf/gDNA_mutant_genes_overlap.pdf", 
    width = 11)
par(mar= c(10,11,6,2),
    tcl= -0.2,
    mgp= c(2,0.5, 0))
dat[, {
  pl <- split(FBgn, cdition)
  vl_upset_plot(pl)
  title(main= class)
  print(class)
}, class]
dev.off()
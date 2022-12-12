setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
genes <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)

dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[type %in% c("nonsynonymous SNV", "stopgain")]
dat[, cdition:= switch(cdition,
                       "PH18_2"= "no ph-KD.1",
                       "PH29_1"= "Constant ph-KD.1",
                       "PH29_2"= "Constant ph-KD.2",
                       "PHD11_1"= "Trans. ph-KD d11.1",
                       "PHD11_2"= "Trans. ph-KD d11.2",
                       "PHD9_1"= "Trans. ph-KD d9.1",
                       "PHD9_2"= "Trans. ph-KD d9.2"), cdition]

pdf("pdf/gDNA_mutant_genes_overlap_with_PH18.pdf", 
    width = 24,
    height = 6)
par(mar= c(12,15,2,0),
    tcl= -0.2,
    mgp= c(2,0.5, 0),
    cxy= 2)
vl_upset_plot(split(dat$id, dat$cdition))
dev.off()

# dat[, .(length(unique(cdition))==6 & !"PH18_2" %in% cdition, symbol), id][(V1)]
# spz3, CG34398 and CG9525
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
dat <- dat[exonic_annotation %in% c("frameshift deletion", "frameshift insertion", "nonsynonymous SNV", "stopgain")]
dat[, cdition:= switch(cdition,
                       "PH18_2"= "No ph-KD.2",
                       "PH29_1"= "Constant ph-KD.1",
                       "PH29_2"= "Constant ph-KD.2",
                       "PHD11_1"= "Trans. ph-KD d11.1",
                       "PHD11_2"= "Trans. ph-KD d11.2",
                       "PHD9_1"= "Trans. ph-KD d9.1",
                       "PHD9_2"= "Trans. ph-KD d9.2"), cdition]
pl <- na.omit(dat[(filter), .(symbol= unlist(symbol)), cdition])

pdf("pdf/Figure_1_mutated_genes_overlap.pdf", 
    width = 20,
    height = 5)
par(mar= c(13,15,2,0),
    tcl= -0.2,
    mgp= c(2,0.5, 0))
vl_upset_plot(split(pl$symbol, pl$cdition), grid_hex = 2.5)
dev.off()

# Find genes mutated in all tumors
sel <- pl[, .("Constant ph-KD.1" %in% cdition
              & "Constant ph-KD.2" %in% cdition
              & "Trans. ph-KD d11.1" %in% cdition
              & "Trans. ph-KD d11.2" %in% cdition
              & "Trans. ph-KD d9.1" %in% cdition
              & "Trans. ph-KD d9.2" %in% cdition
              & !"No ph-KD.2" %in% cdition), symbol][(V1), symbol]
sel <- dat[, .(symbol= unique(unlist(symbol))), 
             .(id, 
               cdition, 
               alt_class,
               allelic_ratio= paste0(TUMOR_ref, ",", TUMOR_alt, " (",
                                     round(TUMOR_alt/(TUMOR_alt+TUMOR_ref)*100), "%)"))][symbol %in% sel]
tab <- dcast(sel, id+symbol~cdition, value.var = "allelic_ratio")
setcolorder(tab, c("id",
                   "symbol",
                   "No ph-KD.2", 
                   "Constant ph-KD.1",
                   "Constant ph-KD.2",
                   "Trans. ph-KD d11.1",
                   "Trans. ph-KD d11.2",
                   "Trans. ph-KD d9.1",
                   "Trans. ph-KD d9.2"))

pdf("pdf/Extended_data_X_mutated_genes_AD.pdf", width = 20)
gridExtra::grid.table(tab)
dev.off()

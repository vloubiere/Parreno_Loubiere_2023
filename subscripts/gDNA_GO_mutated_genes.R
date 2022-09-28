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

pdf("pdf/gDNA_mutant_genes_GO.pdf")
par(mar= c(5,15,5,8),
    tcl= -0.2,
    mgp= c(2,0.5, 0),
    las= 2)
GOs <- dat[, {
  res <- vl_GO_enrich(geneIDs = unique(FBgn), species = "Dm")
  title(main= class)
  res
}, class]
dev.off()

sel <- GOs[padj<0.05 & log2OR>0, 
           .(FBgn= AnnotationDbi::select(org.Dm.eg.db::org.Dm.eg.db,
                                         keys = as.character(GO),
                                         keytype= "GOALL", 
                                         columns= "FLYBASE")$FLYBASE), 
           .(GO, variable)]
sel <- sel[FBgn %in% unlist(dat$FBgn)]
saveRDS(sel, "Rdata/GO_mutated_genes.rds")

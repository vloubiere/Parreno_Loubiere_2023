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
dat <- dat[exonic_annotation %in% c("frameshift deletion", "frameshift insertion", "nonsynonymous SNV", "stopgain")
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
pl <- plot(enr, padj_cutoff= 0.05)
pl <- plot(enr2, padj_cutoff= 0.05)
# par(mar=c(5,20,3,7))
# pl[, {
#   sel <- AnnotationDbi::select(x= org.Dm.eg.db::org.Dm.eg.db,
#                                keys = as.character(variable), # Only GOs from test set are relevant!
#                                keytype= "GOALL", 
#                                columns= "FLYBASE")
#   tab <- dcast(dat[FBgn %in% sel$FLYBASE], symbol~cdition, fun.aggregate = length, drop = F)
#   tab <- tab[order(rowSums(tab[,!"symbol"]), decreasing = T, PH29_1, PH29_2, PHD11_1, PHD11_2, PHD9_1, PHD9_2)]
#   tab <- as.matrix(tab,1)
#   vl_heatmap(tab, cluster_cols= F, display_numbers= T, main= name, col= c("white", "red"))
#   ""
# }, .(variable, name)]
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!(PH18) & type %in% c("nonsynonymous SNV", "stopgain")]
dat <- dat[, .(FBgn= unlist(FBgn), symbol= unlist(symbol)), .(id, class, occurence, cdition)]
dat <- na.omit(dat)
dat[, occurence:= factor(occurence,
                         c("single condition",
                           "shared >=1 conditions"))]
dat[, cdition:= factor(cdition)]
enr <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat$class), species = "Dm")

pdf("pdf/gDNA_mutant_genes_GO.pdf", 4, 4)
par(mar= c(7,31,3,7),
    las= 2,
    cex= 0.5)
pl <- plot(enr)
par(mar=c(5,20,3,7))
pl[, {
  sel <- AnnotationDbi::select(x= org.Dm.eg.db::org.Dm.eg.db,
                               keys = as.character(variable), # Only GOs from test set are relevant!
                               keytype= "GOALL", 
                               columns= "FLYBASE")
  tab <- dcast(dat[FBgn %in% sel$FLYBASE], symbol~cdition, fun.aggregate = length, drop = F)
  tab <- tab[order(rowSums(tab[,!"symbol"]), decreasing = T, PH29_1, PH29_2, PHD11_1, PHD11_2, PHD9_1, PHD9_2)]
  tab <- as.matrix(tab,1)
  vl_heatmap(tab, cluster_cols= F, display_numbers= T, main= name, col= c("white", "red"))
  ""
}, .(variable, name)]
dev.off()

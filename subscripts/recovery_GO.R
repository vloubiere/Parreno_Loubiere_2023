setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]

# Compute GO enrich
GO <- vl_GO_enrich(geneIDs = split(dat$FBgn, dat$recovery), 
                   geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn,
                   species = "Dm",
                   plot= F)

#############################
# PLOT
#############################
pdf("pdf/recovery_GOs.pdf", 4, 5)
par(mar= c(7,31,1,7),
    las= 2,
    cex= 0.5)
plot(GO, 
     padj_cutoff = 0.01,
     cex.balloons= 0.3,
     order= 'log2OR', 
     top_enrich = 20)
par(mar=c(0,0,5,0))
for(patt in c("STAT", "genital"))
{
  GOI <- GO[grepl(patt, name) & padj<0.01, .(name, variable)]
  sel <- AnnotationDbi::select(x= org.Dm.eg.db::org.Dm.eg.db,
                               keys = GOI$variable, # Only GOs from test set are relevant!
                               keytype= "GOALL", 
                               columns= "FLYBASE")
  vl_plot_table(dat[FBgn %in% sel$FLYBASE & recovery=="noRecovery", .(symbol, PRC1_bound, K27me3_bound, K118Ub_bound, recovery)])
  title <- paste0("GOs= ", paste(unique(GOI$name), collapse = ", "))
  title(main= paste0(strwrap(title, width = 50), collapse= "\n"), line= - 3)
}
dev.off()



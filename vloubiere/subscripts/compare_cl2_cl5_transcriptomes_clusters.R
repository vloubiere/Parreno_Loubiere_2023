setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Import data
obj <- readRDS("Rdata/clustering_RNA.rds")
list2env(obj, 
         envir = environment())
# Compute TF FC
dat <- data.table(file= c("db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD11_vs_RNA_WKD.txt",
                          "db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PH29_vs_RNA_W29.txt"),
                  cdition= c("PHD11", "PH29"))
dat <- dat[, fread(file), (dat)]
FC <- dcast(dat, FBgn~cdition, value.var = "log2FoldChange")

sel <- vl_Dmel_motifs_DB_full[gsub("_mot$", "", mot_enrichment_cl2_vs_cl5[padj<0.05, variable]), on= "uniqName_noSpecialChar"]
sel <- merge(sel, 
             FC,
             all.x= T, 
             by= "FBgn",
             sort= F)

# Refind motif analysius cluster 2 vs 5
pdf("pdf/Figures/motifs_cl2_cl5_compaer.pdf", 
    width = 10, 
    height = 18)
par(mar= c(4,10,4,7),
    mfrow= c(1,2))
plot(mot_enrichment_cl2_vs_cl5, 
     padj_cutoff= 0.05)
title(main= "Motif enrich cl 2 vs 5")
vl_heatmap(sel[nrow(sel):1, .(motif_name, PHD11, PH29)],
           rownames = "motif_name",
           cluster_rows=F,
           cluster_cols=F,
           breaks= c(-3,0,3))
dev.off()




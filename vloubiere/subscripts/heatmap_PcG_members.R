setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
dat <- SJ(file= list.files("db/FC_tables/RNA/", "epiCancer_ED_GFP", full.names = T))
dat[, system:= ifelse(grepl("GFP\\+", file), "GFP+", "GFP-")]
dat[, cdition:= gsub(".*system_RNA_(.*).txt$", "\\1", basename(file))]
dat[, cdition:= gsub("RNA_|_RNA", "", cdition)]
dat[, cdition:= factor(cdition, c("PH18_vs_W18", "PHD9_vs_WKD", "PHD11_vs_WKD", "PH29_vs_W29"))]
dat <- dat[, fread(file), (dat)]
gtf <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
gtf <- as.data.table(gtf)
sel <- gtf[gene_symbol %in% c("E(z)", "Caf1-55", "esc", "Su(z)12", "Psc", "Su(z)2", "Pc", "ph-p", "ph-d", "Sce") 
           & type=="gene", .(FBgn= gene_id, symbol= gene_symbol)]
sel[, symbol:= factor(symbol, 
                      c("Pc", "ph-p", "ph-d", "Psc", "Su(z)2", "Sce", "E(z)", "Su(z)12", "esc", "Caf1-55"))]
dat <- dat[sel, on= "FBgn"]

pdf("pdf/Figures/Heatmap_PcG_members_FC.pdf", 
    height = 4,
    width = 6.5)
par(mfrow=c(1,2))
dat[, {
  mat <- dcast(.SD, 
               symbol~cdition, 
               value.var = "log2FoldChange")
  mat <- as.matrix(mat, 1)
  padj <- dcast(.SD, 
                symbol~cdition, 
                value.var = "padj")
  padj <- as.matrix(padj, 1)
  adj <- mat
  adj[!(padj <= 0.01)] <- NA
  vl_heatmap(adj, 
             cluster_rows= F,
             cluster_cols= F, 
             display_numbers= T, 
             main= system, 
             legend_title= "log2FoldChange", 
             breaks= c(-3, 0, 3))
  print(".")
}, system]
dev.off()

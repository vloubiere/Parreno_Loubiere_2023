setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(digest)

dat <- data.table(file= list.files("db/FC_tables/", 
                                   "PHD11_vs_RNA_WKD|PHD9_vs_RNA_WKD|PH29_vs_RNA_W29|EZD11_vs_RNA_WKD|EZD9_vs_RNA_WKD|EZ29_vs_RNA_W29", 
                                   full.names = T))
dat[, cdition:= gsub("RNA_|epiCancer_Parreno_|.txt", "", basename(file)), file]
dat <- dat[, fread(file), (dat)]
dat[padj>=0.01 | is.na(padj), class:= "unaffected"]
dat[padj<0.01 & log2FoldChange>0, class:= "up"]
dat[padj<0.01 & log2FoldChange<0, class:= "down"]
dat[, keep:= ifelse(all(class=="unaffected"), F, T), V1]
dat <- dat[(keep) & class!="unaffected"]

#------------------------#
# COMPARE EZ/PH up/down genes in KD
#------------------------#
sub <- dat[grepl("29.*29$", cdition)]
sub[, cluster:= paste0(cdition, "_", class), .(cdition, class)]
obj <- vl_upset_split(split(sub$V1, sub$cluster))
pl <- obj$elements
names(pl) <- obj$.id

pdf("pdf/genes_overlaps/GO_up_down_genes_EZ_PH_intersections_KD.pdf", 10, 10)
vl_GO_clusters(FBgn_list = pl, 
               N_top = 5, 
               cex = 0.7)
mtext("EZ vs PH GOs KD 29C", line = 0.5)
dev.off()

#------------------------#
# COMPARE EZ/PH up genes KD vs Transient d11
#------------------------#
sub <- dat[class=="up" & grepl("29.*29$|D11.*KD$", cdition)]
obj <- vl_upset_split(split(sub$V1, sub$cdition))
setkeyv(obj, ".id")
pl <- list("EZ+PH KD"= unlist(obj[grepl("EZ29_vs_W29.*PH29_vs_W29", .id), elements]),
           "EZ+PH D11"= unlist(obj[grepl("EZD11_vs_WKD.*PHD11_vs_WKD", .id), elements]))

pdf("pdf/genes_overlaps/GO_up_genes_PH_EZ_KD_vs_D11.pdf", width= 10, height = 10)
vl_GO_clusters(FBgn_list = pl,  
               padj_cutoff = 1e-3, 
               N_top = 10, 
               cex = 0.8)
mtext("UP-REGULATED", line = 0.5)
dev.off()

#------------------------#
# COMPARE EZ/PH DOWN genes KD vs Transient d11
#------------------------#
sub <- dat[class=="down" & grepl("29.*29$|D11.*KD$", cdition)]
obj <- vl_upset_split(split(sub$V1, sub$cdition))
setkeyv(obj, ".id")
pl <- list("EZ+PH KD"= unlist(obj[grepl("EZ29_vs_W29.*PH29_vs_W29", .id), elements]),
           "EZ+PH D11"= unlist(obj[grepl("EZD11_vs_WKD.*PHD11_vs_WKD", .id), elements]))

pdf("pdf/genes_overlaps/GO_down_genes_PH_EZ_KD_vs_D11.pdf", width= 10, height = 10)
vl_GO_clusters(FBgn_list = pl,  
               padj_cutoff = 1e-3, 
               N_top = 10)
mtext("DOWN-REGULATED", line = 0.5)
dev.off()

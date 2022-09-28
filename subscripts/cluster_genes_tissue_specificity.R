setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[is.na(cl), cl:= 0]

# FLYBASE tissue transcriptomes
tsv <- fread("../../genomes/dm6/gene_rpkm_report_fb_2021_02.tsv", fill= T, skip = 5)
colnames(tsv) <- gsub("\\#", "", colnames(tsv))
tsv <- tsv[Parent_library_name == "modENCODE_mRNA-Seq_development"]

dev_hk <- melt(tsv, 
               id.vars = c("FBgn", "GeneSymbol", "RNASource_name"),
               measure.vars = "RPKM_value")
sel <- dev_hk[, sum(value, na.rm= T)>1, FBgn][(V1), FBgn]
dev_hk <- dev_hk[FBgn %in% sel]
setorderv(dev_hk, "value", 1)
dev_hk[, score:= min(which(cumsum(value)>(sum(value)*0.5))), FBgn]
dev_hk[, score:= (score-min(score))/(max(score)-min(score))]
dat[dev_hk, score:= i.score, on="FBgn"]

###########################################
# PLOT
###########################################
pdf("pdf/cluster_PRC1_bound_unbound_tissue_specifity.pdf",
    width= 9,
    height= 6)
par(las= 2,
    mar= c(2,3.5,4,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(3,4))
# Dev hk example
ex <- data.table::copy(vl_genes_set)
ex[dat, score:= i.score, on= "FBgn"]
vl_boxplot(score~GO,
           ex)
# clusters
vl_boxplot(score~PRC1_bound+cl,
           dat,
           col= c("lightgrey", "tomato"),
           compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12), c(13,14)),
           xaxt= "n",
           ylab= "Tissue specificty score",
           at= rep(seq(1, 13, 2), each= 2)+c(0.2, 0.8))
title(main= "Tissue specificity",
      line= 2.75)
abline(h= 0,
       lty= "11")
axis(1, 
     at = seq(1.5, 13.5, 2), 
     labels = paste0("cl ", 0:6))
legend(par("usr")[1],
       par("usr")[4]+strheight("M")*4.5, 
       legend= c("PRC1-", "PRC1+"),
       fill= c("lightgrey", "tomato"),
       bty= "n",
       xpd=T)
dev.off()



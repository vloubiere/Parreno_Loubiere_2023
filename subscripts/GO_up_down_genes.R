setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# Import
dat <- fread("Rdata/final_gene_features_table.txt")
PH18 <- dat[diff_PH18!="unaffected", FBgn]
dat <- melt(dat, id.vars = c("FBgn", "PRC1_bound"), measure.vars = patterns("diff"= "^diff"))
dat <- dat[value!="unaffected" & !(FBgn %in% PH18)]
dat[, cdition:= gsub("diff_", "", variable)]
# GOs
GO <- dat[, {
  vl_GO_enrich(geneIDs = split(FBgn, list(value, ifelse(PRC1_bound, "PRC1+", "PRC1-")), sep= " "),
               species = "Dm", 
               geneUniverse_IDs = fread("Rdata/final_gene_features_table.txt")$FBgn)
}, cdition]

###################################################
# plot
###################################################
pdf("pdf/RNA_GO_up_down_genes_GFP-.pdf", 
    width = 7.5, 
    height = 10)
par(mar= c(8, 27, 4.1, 7),
    las= 2,
    mgp= c(2,0.5,0))
GO[, {
  .c <- .SD
  setattr(.c, "class", c("vl_enr_cl", "data.table", "data.frame"))
  plot(.c,
       padj_cutoff = 0.05,
       top_enrich = 10,
       cex.balloons= 0.45,
       col= c("blue", "red"),
       main= cdition)
  print("done")
}, cdition]
dev.off()

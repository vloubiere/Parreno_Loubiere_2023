# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(ggplot2)
require(ggrepel)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]

#------------------------------#
# Compute and save genes list
#------------------------------#
# Genes list
fwrite(dat[recovery=="Recovery", .(symbol)],
       "db/iCisTarget/Reversible_genes_list.txt",
       col.names = F, 
       row.names = F, 
       sep= "\t",
       quote= F,
       na= NA)
fwrite(dat[recovery=="noRecovery", .(symbol)],
       "db/iCisTarget/Irreversible_genes_list.txt",
       col.names = F, 
       row.names = F, 
       sep= "\t",
       quote= F,
       na= NA)

# Compare genes
dat <- rbind(fread("db/iCisTarget/Irreversible_genes/statistics.tbl", sel= 1:8),
             fread("db/iCisTarget/Reversible_genes/statistics.tbl", sel= 1:8))
NES <- dcast(dat, 
             FeatureID+FeatureDescription+FeatureAnnotations~`#GeneSignatureID`, 
             value.var= "NES")
NES <- NES[FeatureAnnotations!="" & (Irreversible_genes>3 | Reversible_genes>3), 
           .SD[which.max(apply(.SD, 1, max))], FeatureAnnotations, 
           .SDcols= c("Irreversible_genes", "Reversible_genes")]

ggplot(NES, aes(Irreversible_genes, Reversible_genes, label = FeatureAnnotations)) + 
  geom_point(color = "red") + 
  geom_text_repel() +
  geom_abline(intercept = -1, slope = 1)+
  geom_abline(intercept = 1, slope = 1)+
  theme_classic()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- melt(dat[recovery=="noRecovery"], 
            id.vars = c("FBgn", "log2FoldChange_PHD11"),
            measure.vars =  c("PH_PH18_prom", "H3K27me3_PH18_body", "H2AK118Ub_PH18_body"))

pdf("pdf/recovery_PH_levels_vs_FC_D9_D11.pdf", height = 2.75)
par(mfrow= c(1,3))
dat[, {
  PCC <- round(cor.test(log2FoldChange_PHD11, 
                        value)$estimate, 2)
  plot(log2FoldChange_PHD11, 
       value,
       main= paste0("noRecovery genes PCC=", PCC),
       xlab= "gene log2FC PH D11",
       ylab= variable,
       log= "y")
}, variable]
dev.off()
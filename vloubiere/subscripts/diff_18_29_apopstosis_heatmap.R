setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

#-------------------------#
# Apoptosis heatmap
#-------------------------#
apop <- dat[symbol %in% c("hid", "rpr", "Diap1", "Drice", "DCP1", "Decay", "Mmp1", "Dronc", "Dredd", "p53")]
apop <- dcast(apop, 
              symbol~cdition, 
              value.var = "log2FoldChange")
apop <- as.matrix(apop, 1)

pdf("pdf/comparison_ph18_ph29_DOSE/apoptosis_heatmap.pdf", 
    width = 4, 
    height = 7)
vl_heatmap(apop,
           display_numbers = T, 
           legend_title = "log2FoldChange")
dev.off()






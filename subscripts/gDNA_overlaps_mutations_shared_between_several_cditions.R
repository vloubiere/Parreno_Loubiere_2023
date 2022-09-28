setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!(PH18) & occurence=="shared >=1 conditions"]

pdf("pdf/gDNA_overlaps_mutation_shared_several_conditions.pdf", 18, 5)
par(mar= c(10,12,2,2))
dat[, {
  pl <- split(id, cdition)
  vl_upset_plot(pl)
  title(main= class)
}, class]
dev.off()
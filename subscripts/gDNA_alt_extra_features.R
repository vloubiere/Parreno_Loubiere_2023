setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

######################################################################
# Import data
######################################################################
dat <- fread("Rdata/gDNA_final_table.txt")
dat <- dat[!(PH18)]

pdf("pdf/Figures/gDNA_conditions_overlaps.pdf", 
    15, 
    4.5)
par(cex= 0.6)
dat[, {
  ids <- split(id, cdition)
  vl_upset_plot(ids)
  title(main= paste(class, alt_class))
}, .(class, alt_class)]
dev.off()

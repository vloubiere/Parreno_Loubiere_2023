setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(vlfunctions)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat[cdition=="PH18_2", cdition:= "PH18"]
dat[, cdition:= gsub("_", "_rep", cdition)]
dat[, cdition:= factor(cdition, 
                       c("PH18",
                         "PH29_rep1",
                         "PH29_rep2",
                         "PHD9_rep1",
                         "PHD9_rep2",
                         "PHD11_rep1",
                         "PHD11_rep2"))]
setorderv(dat, "cdition")

pdf("pdf/gDNA_conditions_overlaps.pdf", 
    30, 
    4.5)
par(cex= 0.6,
    mar= c(12,15,2,2))
dat[, {
  pl <- split(id, cdition)
  vl_upset_plot(pl)
  title(main=  alt_class)
}, alt_class]
dev.off()
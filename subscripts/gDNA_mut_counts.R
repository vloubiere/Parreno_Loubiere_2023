setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(vlfunctions)

######################################################################
# Import data
######################################################################
dat <- readRDS("Rdata/gDNA_final_table.rds")
dat <- dat[!(PH18)]
dat[, cdition:= gsub("_", "_rep", cdition)]
dat[, cdition:= factor(cdition, 
                       c("PH29_rep1",
                         "PH29_rep2",
                         "PHD9_rep1",
                         "PHD9_rep2",
                         "PHD11_rep1",
                         "PHD11_rep2"))]
setorderv(dat, "cdition")

pdf("pdf/gDNA_mutations_counts.pdf", 4, 4)
dat[, {
  vl_plot_table(.SD[, .(Count= .N), cdition])
  title(main= paste0(alt_class, "\n(not found in PH18)"))
}, alt_class]
dat[, {
  total <- length(unique(id))
  counts <- .SD[, {
    .c <- length(unique(id))
    .(Count= paste0(.c, " (", round(.c/total*100, 1), "%)"))
  }, occurence]
  vl_plot_table(rbind(counts,
                      data.table(occurence= "Total", Count= total)))
  title(main= paste0(alt_class, "\n(not found in PH18)"))
}, alt_class]
dev.off()



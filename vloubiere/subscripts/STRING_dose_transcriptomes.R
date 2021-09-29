setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")


#----------------------#
# Download stringdb
#----------------------#
dat <- readRDS("Rdata/clustering_dose_transcriptomes.rds")
agg <- dat[, .(value= mean(value)), .(rcl, symbol)]

#----------------------#
# Plot networks
#----------------------#
pdf("pdf/clustering/clustering_dose_STRING.pdf")
agg[, {
  interactions <- vl_STRING_interaction(symbols = symbol, 
                                        size = abs(value),
                                        col = ifelse(value<0, "cornflowerblue", "tomato"),
                                        score_cutoff = switch(rcl,
                                                              "1"= 400,
                                                              "2"= 400,
                                                              "3"= 800,
                                                              "4"= 400,
                                                              "5"= 400))
  lab.cex <- interactions$vertices$size
  lab.cex <- lab.cex/max(lab.cex)
  vl_STRING_network(obj = interactions, 
                    cex.vertices = 3, 
                    label.cex = lab.cex*0.3)
  mtext(paste("Cluster #", rcl), side= 1)
  print(paste0(rcl, " DONE!"))
}, rcl]
dev.off()

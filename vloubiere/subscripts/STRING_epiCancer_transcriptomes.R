setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")


#----------------------#
# Download stringdb
#----------------------#
dat <- readRDS("Rdata/clustering_epiCancer_transcriptomes.rds")
agg <- dat[, .(value= mean(value)), .(rcl, symbol)]

#----------------------#
# Plot networks
#----------------------#
pdf("pdf/clustering/clustering_epiCancer_STRING.pdf")
agg[, {
  interactions <- vl_STRING_interaction(symbols = symbol, 
                                        size = abs(value),
                                        col = ifelse(value<0, "cornflowerblue", "tomato"),
                                        score_cutoff = switch(rcl,
                                                              "1"= 700,
                                                              "2"= 100,
                                                              "3"= 400,
                                                              "4"= 200,
                                                              "5"= 200,
                                                              "6"= 200,
                                                              "7"= 200,
                                                              "8"= 200), 
                                        cex.label = abs(value))
  vl_STRING_network(obj = interactions, 
                    cex.vertices = 5, 
                    cex.vertices.labels = 0.6)
  mtext(paste("Cluster #", rcl), side= 1)
  print(paste0(rcl, " DONE!"))
}, rcl]
dev.off()

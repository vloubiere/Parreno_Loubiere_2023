setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- melt(dat, 
            id.vars = c("FBgn", "recovery", "cl", "seqnames", "start", "end", "strand"),
            measure.vars = c("FPKM_PH18", "FPKM_PH29", "FPKM_PHD9", "FPKM_PHD11",
                             "H3K27Ac_PH18_prom", "H3K27Ac_PH29_prom", "H3K27Ac_PHD9_prom", "H3K27Ac_PHD11_prom"))
dat[, c("variable", "cdition"):= tstrsplit(variable, "_", keep= 1:2)]
dat <- dcast(dat, FBgn+seqnames+start+end+strand+recovery+cl+cdition~variable, value.var = "value")

# control set
ctlSet <- dat[!is.na(recovery), {
  ctls <- dat[is.na(cl)][.BY, !"recovery", on= "cdition"]
  idx <- seq(.N)
  res <- .SD[, {
    c_FPKM <- FPKM
    .c <- ctls[order(abs(FPKM-c_FPKM))][1]
    ctls <- ctls[FBgn!=.c$FBgn]
    .c
  }, idx]
  res$cdition <- res$idx <- NULL
  res
}, .(cdition, recovery)]
ctlSet[, recovery:= paste0(gsub("ecovery", "", recovery), "_control")]

# Format and resize
dat <- rbind(dat[!is.na(recovery)], 
             ctlSet)
dat[, track:= list.files("db/bw/cutnrun/", paste0("H3K36me3_", cdition, "_merge"), full.names = T), cdition]
dat <- vl_resizeBed(dat, center = "start", upstream = 0, downstream = 0)
dat[, recovery:= factor(recovery)]

# PLOT
Cc <- c("rosybrown1", "rosybrown3", "palegreen1", "palegreen3")

pdf("pdf/recovery_H3K27Ac_activity_matched_genes.pdf",
    width = 6, 
    height = 9)
layout(matrix(1:12, 4, 3, byrow = T), 
       widths = c(0.35,0.35,0.9))
par(mar= c(5,3,2,1),
    mgp= c(1.5,0.5,0),
    tcl= -0.2,
    las= 1)
dat[, {
  vl_boxplot(FPKM~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2), c(3,4)),
             main= cdition,
             col= Cc,
             ylab= "FPKM")
  vl_boxplot(H3K27Ac~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2), c(3,4)),
             main= cdition,
             col= Cc,
             ylab= "H3K27Ac enrichment")
  vl_bw_average_track(bed = copy(.SD), 
                      tracks = track, 
                      upstream = 2500,
                      downstream = 10000,
                      set_IDs = recovery, 
                      stranded = T, 
                      legend_pos = "topright",
                      center_label = "TSS", 
                      col= Cc)
  abline(v= c(-750, 250), lty= 2)
  print("")
}, .(cdition, track)]
dev.off()
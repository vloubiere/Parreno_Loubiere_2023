setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- melt(dat, 
            id.vars = c("FBgn", "recovery", "cl", "seqnames", "start", "end", "strand"),
            measure.vars = c("FPKM_PH18", "FPKM_PH29", "FPKM_PHD9", "FPKM_PHD11",
                             "H3K36me3_PH18_TTS", "H3K36me3_PH29_TTS", "H3K36me3_PHD9_TTS", "H3K36me3_PHD11_TTS"))
dat[, c("variable", "cdition"):= tstrsplit(variable, "_", keep= 1:2)]
dat <- dcast(dat, FBgn+seqnames+start+end+strand+recovery+cl+cdition~variable, value.var = "value")
dat[, gene_width:= end-start+1]

# control set
ctlSet <- dat[recovery=="noRecovery", {
  ctls <- dat[is.na(cl)][.BY, on= "cdition"]
  idx <- seq(.N)
  res <- .SD[, {
    c_FPKM <- FPKM   
    c_width <- gene_width
    size_check <- ctls[, between(gene_width, c_width*0.95, c_width*1.05)]
    .c <- ctls[size_check][order(abs(FPKM-c_FPKM))][1]
    ctls <- ctls[FBgn!=.c$FBgn]
    .c
  }, idx]
  res$cdition <- res$idx <- NULL
  res
}, cdition]
ctlSet[, recovery:= "control"]

# Format and resize
dat <- rbind(dat[recovery=="noRecovery"], 
             ctlSet)
dat[, track:= list.files("db/bw/cutnrun/", paste0("H3K36me3_", cdition, "_merge"), full.names = T), cdition]
dat <- vl_resizeBed(dat, center = "end", upstream = 0, downstream = 0)

# PLOT
Cc <- c("lightgrey", "rosybrown1")

pdf("pdf/recovery_H3K36me3_noRecovery_genes.pdf", width = 6, height = 9)
layout(matrix(1:16, 4, 4, byrow = T), 
       widths = c(0.35,0.35,0.35,1))
par(mar= c(5,3,2,1),
    mgp= c(1.5,0.5,0),
    tcl= -0.2,
    las= 1)
dat[, {
  vl_boxplot(FPKM~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2)),
             main= cdition,
             col= Cc,
             ylab= "FPKM")
  vl_boxplot(log10(gene_width+1)~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2)),
             main= cdition,
             col= Cc,
             ylab= "log10(Gene length+1)")
  vl_boxplot(H3K36me3~recovery,
             tilt.names= T,
             compute_pval= list(c(1,2)),
             main= cdition,
             col= Cc,
             ylab= "H3K36me3 enrichment")
  vl_bw_average_track(bed = copy(.SD), 
                      tracks = track, 
                      upstream = 10000,
                      downstream = 10000,
                      set_IDs = recovery, 
                      stranded = T, 
                      center_label = "TTS", 
                      col= Cc)
  abline(v= c(-2500, 1000), lty= 2)
  print("")
}, .(cdition, track)]
dev.off()
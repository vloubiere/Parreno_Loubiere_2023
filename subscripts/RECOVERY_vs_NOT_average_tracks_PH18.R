setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(recovery)]

# Melt
.m <- melt(dat, id.vars = "recovery", measure.vars = patterns("_body$|_prom$"))
.m <- .m[grepl("PH18|ED", variable)]
.m[, variable:= gsub("_body$|_prom$", "", variable)]
.m[, track:= list.files("db/bw/cutnrun/", 
                        paste0(variable, "_merge"), 
                        full.names = T), variable]
.m[is.na(track), track := list.files("db/bw/SA_2020/", 
                                     paste0(gsub("_body$|_prom$", "", variable), "_merge"), 
                                     full.names = T), variable]

# Average tracks window
gtf <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)[type=="gene"]
TSS <- vl_resizeBed(gtf, "start", 0, 0)
TSS[dat, recovery:= i.recovery, on="gene_id==FBgn"]
TSS <- na.omit(TSS[, .(seqnames, start, end, strand, recovery)])
TSS[, recovery:= factor(recovery, c("Recovery", "noRecovery"))]

#-----------------------------------------#
# Plot
#-----------------------------------------#
Cc <- c("palegreen3", "rosybrown1")
pdf("pdf/Figures/PH18_CUTNRUN_enrich_revert_vs_not.pdf",
    height = 5, 
    width = 30)
layout(matrix(1:(6*4), nrow=2, byrow = T), 
       widths= rep(c(1,0.5), 12))
par(mar= c(5,4,2,1),
    las= 1,
    mgp= c(2.5,0.5,0),
    tcl= -0.2)
.m[, {
  # Plot average track for each class
  .q <- vl_bw_average_track(bed= TSS, 
                            tracks= track,
                            set_IDs = TSS$recovery,
                            upstream = 2500,
                            downstream = 10000,
                            stranded = T,
                            center_label = "TSS", 
                            legend.cex = 0.6, 
                            legend= F,
                            col = Cc)
  
  # Legend
  legend("topright",
         legend = c(paste0("Recovery (", sum(recovery=="Recovery"), " genes)"),
                    paste0("No recovery (", sum(recovery=="noRecovery"), " genes)")),
         fill= Cc,
         bty= "n")
  title(main= variable)
  
  # Boxplot
  vl_boxplot(value~recovery,
             col= Cc,
             ylab= "Enrichment", 
             compute_pval= list(c(1,2)),
             tilt.names= T)
  print(".")
}, .(variable, track)]
dev.off()

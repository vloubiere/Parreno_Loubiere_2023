setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/final_gene_features_table.txt")
K27me3 <- vl_importBed("db/peaks/cutnrun/H3K27me3_PH18_confident_peaks.broadPeak")
K27me3[, HOX:= grepl("_1041$|_1040$|_901$|_907$|_909$", name)]
K27me3[, gp:= vl_collapseBed(K27me3, 10000, return_idx_only = T)]
K27me3 <- K27me3[, .(start= min(start), end= max(end), HOX= any(HOX)), .(seqnames, gp)]
TSS <- vl_resizeBed(dat, "start", 0, 0)
cl <- vl_closestBed(TSS, K27me3)
dat$domain_size_all <- cl[, end.b-start.b]/1000
cl <- vl_closestBed(TSS, K27me3[!(HOX)]) # kikck out BXC ANTP
dat$domain_size_noBX <- cl[, end.b-start.b]/1000

###########################################
# PLOT
###########################################
pdf("pdf/recovery_size_K27_domains.pdf", 
    width = 5,
    height = 5)
par(las= 1,
    mar= c(8,10,8,10),
    tcl= -0.2,
    mgp= c(2,0.5,0))
vl_boxplot(domain_size_all~recovery,
           dat[!is.na(recovery)],
           col=  c("rosybrown1", "palegreen3"),
           ylab= "Length K27me3 domain in kb",
           compute_pval= list(c(1,2)),
           notch= T,
           tilt.names= T,
           main= "all domains")
vl_boxplot(domain_size_noBX~recovery,
           dat[!is.na(recovery)],
           col=  c("rosybrown1", "palegreen3"),
           ylab= "Length K27me3 domain in kb",
           compute_pval= list(c(1,2)),
           notch= T,
           tilt.names= T,
           main= "No BX-C ANTP")
dev.off()
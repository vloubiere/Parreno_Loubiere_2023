setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(vlfunctions)
require(GenomicRanges)

#---------------------------------------------------------------#
# H3K27Ac
#---------------------------------------------------------------#
dat <- as.data.table(readRDS("Rdata/K27Ac_cutnrun_clustering.rds"))
dat <- dat[, max(abs(value)), .(rn, rn_cl)][, .SD[order(V1, decreasing = T)][1:5], rn_cl]
setorderv(dat, "rn_cl")

tracks <- data.table(file= c("db/bw/SA_2020/PH_ED_merge.bw",
                             "db/bw/SA_2020/EZ_CNSID_BL_merge.bw",
                             "db/bw/SA_2020/SUZ12_ED_merge.bw",
                             normalizePath(list.files("db/bw/cutnrun_merge/", "H3K27Ac", full.names = T))))

pdf("pdf/cutnrun/screenshots_K27Ac_diff_cluters.pdf", 
    width = 10, 
    height = 6)
dat[, {
  .c <- GRanges(rn)
  .regions <- resize(.c, width = width(.c)+5000, fix = "center")
  par(mar= c(0.5, 10.1, 2.1, 0.5))
  vl_screenshot(.regions,
                max = c(4000,4000,2000,2500,2500,2500,2500),
                highlight_regions = .c,
                tracks$file, 
                genome = "dm6", 
                n_genes = 1)
  mtext(paste0("H3K27Ac cluster= ", rn_cl), cex= 2)
}, rn_cl]
dev.off()

#---------------------------------------------------------------#
# H3K27me3
#---------------------------------------------------------------#
dat <- as.data.table(readRDS("Rdata/K27me3_cutnrun_clustering.rds"))
dat <- dat[, max(abs(value)), .(rn, rn_cl)][, .SD[order(V1, decreasing = T)][1:5], rn_cl]
dat <- na.omit(dat)
setorderv(dat, "rn_cl")

tracks <- data.table(file= c("db/bw/SA_2020/PH_ED_merge.bw",
                             "db/bw/SA_2020/EZ_CNSID_BL_merge.bw",
                             "db/bw/SA_2020/SUZ12_ED_merge.bw",
                             normalizePath(list.files("db/bw/cutnrun_merge/", "H3K27me3", full.names = T))))

pdf("pdf/cutnrun/screenshots_K27me3_diff_cluters.pdf", 
    width = 10, 
    height = 6)
dat[, {
  .c <- GRanges(rn)
  .regions <- resize(.c, width = width(.c)+5000, fix = "center")
  par(mar= c(0.5, 10.1, 2.1, 0.5))
  vl_screenshot(.regions,
                max = c(4000,4000,2000,2500,2500,2500,2500),
                highlight_regions = .c,
                tracks$file, 
                genome = "dm6", 
                n_genes = 1)
  mtext(paste0("H3K27me3 cluster= ", rn_cl), cex= 2)
}, rn_cl]
dev.off()

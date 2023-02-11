# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import counts from featureCounts file
dat <- readRDS("db/counts/RNA/epiCancer_GFP_GFP_counts.rds")
counts <- as.data.table(dat[[1]], keep.rownames = "rn")
counts <- melt(counts, id.vars = "rn")
counts <- counts[rn %in% c("EGFP", "mRFP1") & grepl("GFP.PH", variable)]
tot <- as.data.table(t(dat$stat[12,-1]), keep.rownames = T)
setnames(tot, "12", "total")
counts[tot, total:= i.total, on= "variable==rn"]
counts[, cdition:= tstrsplit(variable, "[.]", keep= 3)]
counts[, cdition:= switch(cdition, 
                          "PH18"= "No ph-KD",
                          "PH29"= "Constant ph-KD",
                          "PHD9"= "Day9 ",
                          "PHD11"= "Day 11"), cdition]
counts[, cdition:= factor(cdition, c("No ph-KD","Constant ph-KD","Day9 ","Day 11"))]
# libsize norm
counts[, norm_counts:= value/total*1e6]
# Select GFP form
counts <- counts[, .(mean= mean(value), sd= sd(value), value= list(value)), .(rn, cdition)]
counts[, col:= fcase(rn=="mRFP1", "tomato", default= "limegreen")]

pdf("pdf/RNA_GFP_norm_counts.pdf", 
    width = 1.9,
    height= 8)
par(mfrow= c(3,1),
    mar= c(10,6,4,3),
    mgp= c(3,0.5,0),
    tcl= -0.2)
counts[, {
  bar <- barplot(mean,
                 las= 2,
                 names.arg = cdition,
                 col= col, 
                 ylab= paste0(rn,  " normalized counts"),
                 ylim= c(0, c(8e3, 1e5, 2e4)[.GRP]),
                 xaxt= "n")
  vl_tilt_xaxis(bar, labels= levels(cdition))
  points(rep(bar, each= 3),
         unlist(value),
         xpd= T,
         pch= 16,
         col= "darkgrey",
         cex= 0.7)
  arrows(bar, mean, bar, mean+sd, length = 0.02, angle = 90)
  arrows(bar, mean, bar, mean-sd, length = 0.02, angle = 90)
  .SD
}, rn]
dev.off()

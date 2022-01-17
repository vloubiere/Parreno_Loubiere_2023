setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(ggalluvial)

#-------------------#
dat <- rbindlist(readRDS("Rdata/RNA_tables_object.rds")$FC, idcol = T)
dat[padj<0.05 & log2FoldChange>1, class:="up"]
dat[padj<0.05 & log2FoldChange<(-1), class:="down"]
dat[is.na(class), class:="unaffected"]
dat <- dat[!(cdition=="PHD9" & .id=="allograft")]
dat[, cdition:= factor(cdition, levels= c("PH18",
                                          "PHD11", 
                                          "PHD11_T2", 
                                          "PHD11_T3",
                                          "PHD9",
                                          "PH29"))]
                                  

pdf("pdf/RNA_timecourse/alluvial_plot_timecourse.pdf", width = 14, height = 4)
par(mfrow=c(1,2), 
    mar= rep(1,4),
    cex= 0.7)
dat[, {
  .c <- dcast(.SD, FBgn~cdition, value.var = "class")[, -1]
  .c <- .c[apply(.c, 1, function(x) any(x!="unaffected"))]
  vl_alluvial_plot(.c,
                   class_levels = c("down", "unaffected", "up"))
}, .id]
dev.off()

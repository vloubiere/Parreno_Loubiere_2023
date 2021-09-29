setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(ggalluvial)

#-------------------#
dat <- readRDS("Rdata/final_FC_table.rds")
dat[padj<0.01 & log2FoldChange>0, class1:="up"]
dat[padj<0.01 & log2FoldChange<0, class1:="down"]
dat[is.na(class1), class1:="unaffected"]
dat[padj<0.01 & log2FoldChange>1, class2:="up"]
dat[padj<0.01 & log2FoldChange<(-1), class2:="down"]
dat[is.na(class2), class2:="unaffected"]

pdf("pdf/alluvial_plot_timecourse.pdf", width = 14, height = 4)
par(mfrow=c(1,2), 
    mar= rep(1,4))
# PH dose low cutoff
dose <- dat[grepl("dose$", cdition), .(any(class1!="unaffected"), cdition, class1), FBgn][(V1)]
vl_alluvial_plot(dcast(dose, FBgn~cdition, value.var = "class1")[, -1],
                 class_levels = c("down", "unaffected", "up"))
# Dose high cutoff
dose <- dat[grepl("dose$", cdition), .(any(class2!="unaffected"), cdition, class2), FBgn][(V1)]
vl_alluvial_plot(dcast(dose, FBgn~cdition, value.var = "class2")[, -1],
                 class_levels = c("down", "unaffected", "up"))

# PH TS low cutoff
PH_TS <- dat[grepl("PH.*TS$", cdition), .(any(class1!="unaffected"), cdition, class1), FBgn][(V1)]
vl_alluvial_plot(dcast(PH_TS, FBgn~cdition, value.var = "class1")[, -1], 
                 class_levels = c("down", "unaffected", "up"))

# PH TS high cutoff
PH_TS <- dat[grepl("PH.*TS$", cdition), .(any(class2!="unaffected"), cdition, class2), FBgn][(V1)]
vl_alluvial_plot(dcast(PH_TS, FBgn~cdition, value.var = "class2")[, -1], 
                 class_levels = c("down", "unaffected", "up"))


# EZ TS low cutoff
EZ_TS <- dat[grepl("EZ.*TS$", cdition), .(any(class1!="unaffected"), cdition, class1), FBgn][(V1)]
vl_alluvial_plot(dcast(EZ_TS, FBgn~cdition, value.var = "class1")[, -1], 
                 class_levels = c("down", "unaffected", "up"))

# EZ TS high cutoff
EZ_TS <- dat[grepl("EZ.*TS$", cdition), .(any(class2!="unaffected"), cdition, class2), FBgn][(V1)]
vl_alluvial_plot(dcast(EZ_TS, FBgn~cdition, value.var = "class2")[, -1], 
                 class_levels = c("down", "unaffected", "up"))

dev.off()


setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/final_gene_features_table.txt")
dat[, cl:= as.character(cl)]

.m <- melt(dat[!is.na(recovery)], id.vars = "recovery", measure.vars = patterns("log2FoldChange"))
.m[, variable:= gsub("^log2FoldChange_", "", variable)]

###########################################
# PLOT
###########################################
pdf("pdf/recovery_FC_conditions.pdf", 
    width = 5,
    height = 5)
par(las= 1,
    mar= c(7,10,3,7),
    tcl= -0.2,
    mgp= c(2,0.5,0))
test <- vl_boxplot(value~recovery+variable, 
           .m, ylab= "log2FoldChange",
           col= c("rosybrown1", "palegreen3"),
           compute_pval= list(c(1,2), c(3,4), c(5,6), c(7,8)),
           xaxt= "n",
           tilt.names= T,
           at= rep(seq(1,7,2), each= 2)+c(0.2,0.8))
axis(1, 
     at= c(1.5, 3.5, 5.5, 7.5),
     c("PH18", "PH29", "PHD11", "PHD9"),
     las= 2)
legend(par("usr")[2],
       par("usr")[4],
       fill= c("palegreen3", "rosybrown1"),
       legend= c("Recovery", "No recovery"),
       bty= "n",
       xpd= T)
dev.off()
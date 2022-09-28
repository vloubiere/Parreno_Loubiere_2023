setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import and compute features
#############################
# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat <- dat[!is.na(cl)]
dat <- dat[!is.na(recovery)]
dat[, col:= ifelse(recovery=="Recovery", "palegreen3", "rosybrown1")]

# Compute gene network
cols <- c("log2FoldChange_PH29", "log2FoldChange_PHD9", "log2FoldChange_PHD11")
sizes <- apply(dat[, ..cols], 1, function(x) max(abs(x), na.rm=T))
sizes[sizes>quantile(sizes, 0.90)] <- quantile(sizes, 0.90)
sizes <- log2(sizes+2)
net <- vl_STRING_interaction(symbols = dat$symbol,
                             species = "Dm",
                             col= dat$col,
                             size = sizes*2.5,
                             cex.label = sizes/4,
                             plot= F)

#############################
# PLOT
#############################
pdf("pdf/recovery_STRING_network.pdf", 3, 3)
par(mar= c(1,1,1,1),
    las= 1,
    cex= 0.5)
set.seed(1)
plot(net)
legend("topleft",
       fill= c("palegreen3", "rosybrown1"),
       legend= c("Recovery", "No recovery"),
       bty= "n",
       border= NA,
       cex= 1.5)
dev.off()

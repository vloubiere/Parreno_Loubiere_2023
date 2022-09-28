setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import counts from featureCounts file!
meta <- fread("Rdata/processed_metadata_RNA.txt")
counts <- meta[, {
    .c <- readRDS(GFP_counts)
    res <- as.data.table(.c[[1]], keep.rownames = "rn")
    res <- melt(res, id.vars = "rn")
    tot <- as.matrix(.c[[4]][,-1], 1)
    res[, tot:= sum(tot[, variable]), variable]
    res
}, .(system, GFP_counts)]
counts[, cdition:= tstrsplit(variable, "[.]", keep= 3)]
# libsize norm
counts[, norm_counts:= value/tot*1e6]
# Select GFP form
counts <- counts[(system=="GFP" & rn!="GFP") | (system=="noGFP" & rn=="GFP")]
counts <- counts[, .(mean= mean(value), value= list(value)), .(rn, system, cdition)]
counts[, col:= fcase(rn=="mRFP1", "tomato", default= "limegreen")]

pdf("pdf/RNA_GFP_norm_counts.pdf", 
    width = 2.1,
    height= 9)
par(mfrow= c(3,1),
    mar= c(10,6,4,3),
    mgp= c(3,0.5,0),
    tcl= -0.2)
counts[, {
  bar <- barplot(mean,
                 las= 2,
                 names.arg = cdition,
                 col= col, 
                 main = system,
                 ylab= paste0(rn,  " normalized counts"),
                 ylim= c(0, c(8e3, 1e5, 2e4)[.GRP]))
  points(rep(bar, each= 3),
         unlist(value),
         xpd= T,
         pch= 16,
         col= "darkgrey",
         cex= 0.7)
}, .(rn, system)]
dev.off()

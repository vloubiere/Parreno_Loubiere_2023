setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

counts <- fread("db/counts/GFP_RFP/GFP_RFP_counts_epiCancer_RNA.txt")
counts <- counts[, .N, .(cdition, dds_file, DESeq2_object, bam, V1)]
sizeFactors <- counts[, as.data.table(readRDS(dds_file)$sizeFactor, keep.rownames = T), dds_file]
sizeFactors[, V1:= gsub(".bam", "_GFP.bam", V1)]
counts[, bam_basename:= basename(bam)]
counts[sizeFactors, sizeFactor:= i.V2, on= c("bam_basename==V1", "dds_file")]
counts <- na.omit(counts[, .(norm_counts= .(N/sizeFactor),
                             mean= mean(N/sizeFactor)), .(DESeq2_object, cdition, V1)])
counts[, genotype:= fcase(grepl("GFP\\+", DESeq2_object), "GFP+ genotype",
                          grepl("GFP-", DESeq2_object), "GFP- genotype"), DESeq2_object]
counts[, col:= fcase(grepl("GFP", V1), "limegreen",
                     grepl("RFP", V1), "tomato"), DESeq2_object]

pdf("pdf/Figures/GFP_norm_counts.pdf", 
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
                 main = genotype,
                 ylab= paste0(V1,  " normalized counts"),
                 ylim= c(0, c(1e4, 1e5, 2.5e4)[.GRP]))
  points(rep(bar, each= 3),
         unlist(norm_counts),
         xpd= T,
         pch= 16,
         col= "darkgrey",
         cex= 0.7)
}, .(DESeq2_object, genotype, V1, col)]
dev.off()
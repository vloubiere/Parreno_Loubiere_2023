setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
# require(vlfunctions)

#-------------------#
dat <- readRDS("Rdata/final_FC_table.rds")
dat[padj<0.01 & log2FoldChange>0, class1:= "up"]
dat[padj<0.01 & log2FoldChange<0, class1:= "down"]
dat[padj<0.01 & log2FoldChange>1, class2:= "up"]
dat[padj<0.01 & log2FoldChange<(-1), class2:= "down"]
.l1 <- split(dat[!is.na(class1), FBgn], dat[!is.na(class1), cdition])
.l2 <- split(dat[!is.na(class2), FBgn], dat[!is.na(class2), cdition])

# Chart
pdf("pdf/genes_overlaps/upsetplot_overlap_up_down_genes.pdf", 
    width = 20)
par(mfrow= c(1,2))

# Dose only
vl_upset_plot(.l1[grep("dose$", names(.l1))], 
              intersection_cutoff = 50)
vl_upset_plot(.l2[grep("dose$", names(.l1))],
              intersection_cutoff = 50)

# PH TS
vl_upset_plot(.l1[grep("PH.*TS$", names(.l1))])
vl_upset_plot(.l2[grep("PH.*TS$", names(.l1))])

# EZ TS
vl_upset_plot(.l1[grep("EZ.*TS$", names(.l1))])
vl_upset_plot(.l2[grep("EZ.*TS$", names(.l1))])

# ALL TS
vl_upset_plot(.l1[grep("TS$", names(.l1))], 
              intersection_cutoff = 100)
vl_upset_plot(.l2[grep("TS$", names(.l1))], 
              intersection_cutoff = 50)

# All
vl_upset_plot(.l1, 
              intersection_cutoff = 100)
vl_upset_plot(.l2, 
              intersection_cutoff = 50)

dev.off()
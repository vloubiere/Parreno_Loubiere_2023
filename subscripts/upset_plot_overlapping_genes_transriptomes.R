setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import data ----
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- melt(dat,
            "FBgn",
            c("diff_PH18",
              "diff_PH29",
              "diff_PHD9",
              "diff_PHD11"))
dat <- dat[value!="unaffected"]
dat[, name:= switch(as.character(variable), 
                    "diff_PH18"= "No ph-KD",
                    "diff_PH29"= "Constant ph-KD",
                    "diff_PHD9"= "Transient ph-KD d9",
                    "diff_PHD11"= "Transient ph-KD d11",), variable]

# Plot ----
pdf("pdf/RNA_upsetplot_overlapping_genes.pdf", 5, 5)
vl_par(mai= c(2.5, 2.75, .2, .2),
       mgp= c(1, .5, 0))
vl_upset_plot(split(dat$FBgn, dat[, .(name, value)]), 
              intersection.cutoff = 10,
              cex.grid = .8,
              grid.hex = .7,
              cex.inter = 5/12)
dev.off()

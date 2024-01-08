setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import clustering data and extra features ----
dat <- readRDS("Rdata/clustering_RNA_GOs.rds") # Extra
dat <- dat[set_hit>2 & name %in% dat[padj<0.01, name]]
setorderv(dat, "log2OR", -1)
dat[, rank:= rowid(cl)]
dat[, cl:= factor(cl,
                  c("Reversible PcG+",
                    "Irreversible PcG+",
                    "Transient-specific PcG+",
                    "Down 1 PcG+",
                    "Down 2 PcG+",
                    "Down 3 PcG+",
                    "Reversible PcG-",
                    "Irreversible PcG-",
                    "Transient-specific PcG-",
                    "Down 1 PcG-",
                    "Down 2 PcG-",
                    "Down 3 PcG-"))]

# select ----
sel <- c("DNA replication",
         "DNA recombination",
         "DNA repair",
         "paracrine signaling",
         "cytokine activity",
         "anterior/posterior pattern specification, imaginal disc",
         "receptor signaling pathway via JAK-STAT",
         "regulation of cell population proliferation",
         "DNA-binding transcription factor activity",
         "cell fate commitment",
         "structural constituent of cuticle",
         "eye-antennal disc morphogenesis",
         "head segmentation",
         "cell fate specification",
         "transmembrane transporter activity",
         "cellular response to calcium ion",
         "segmentation",
         "eye photoreceptor cell differentiation",
         "eye development",
         "neuron differentiation",
         "posterior head segmentation",
         "anterior/posterior pattern specification",
         "cellular response to calcium ion",
         "dopamine secretion",
         "synaptic signaling",
         "compound eye cone cell differentiation",
         "eye photoreceptor cell fate commitment")
if(!all(sel %in% dat$name))
  stop("Some selecting GOs missing!")

# Plot ----
pdf("pdf/Figure_2_clusters_selected_GOs.pdf", 
    width= 8.3,
    height= 11.7)
par(las= 2,
    mar= c(21,26,25,9),
    lwd= 0.5,
    cex.axis= 0.6)
plot(dat[name %in% sel],
     order= "log2OR",
     padj.cutoff = 0.05,
     cex.balloons= 0.2)
title("GOs enrichment per cluster +/- PcG")
dev.off()
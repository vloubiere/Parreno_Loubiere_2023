# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

#############################
# Import clustering data and extra features
#############################
obj <- readRDS("Rdata/clustering_RNA_features.rds") # Extra
GO_PcG <- obj$GO_PcG

#################################
# PLOT
#################################
pdf("pdf/Figure_2_clusters_all_GOs.pdf", 
    width= 8.3,
    height= 11.7)
par(las= 2,
    mar= c(4,30,3,7),
    lwd= 0.5,
    cex.axis= 0.3)
sel <- plot(GO_PcG,
            padj_cutoff = 0.05,
            top_enrich= 50,
            cex.balloons= 0.1)
title("GOs enrichment per cluster +/- PcG")
dev.off()


sel <- sel[name %in% c("DNA replication",
                       "cellular response to DNA damage stimulus",
                       "DNA repair",
                       "MCM complex",
                       "mitotic cell cycle",
                       "cellular developmental process",
                       "tissue development",
                       "sequence−specific DNA binding",
                       "imaginal disc development",
                       "cytokine activity",
                       "anterior/posterior pattern specification, imaginal disc",
                       "receptor signaling pathway via JAK-STAT",
                       "segment specification",
                       "structural constituent of cuticle",
                       "extracellular matrix",
                       "contractile fiber",
                       "myofilament",
                       "sensory organ development",
                       "eye development",
                       "eye−antennal disc development",
                       "eye photoreceptor cell differentiation",
                       "regulation of cellular biosynthetic process",
                       "cell projection",
                       "neuron projection",
                       "axon",
                       "cell junction",
                       "regulation of trans−synaptic signaling",
                       "synapse organization",
                       "cell junction organization",
                       "eye photoreceptor cell fate commitment")]

pdf("pdf/Figure_2_clusters_selected_GOs.pdf", 
    width= 8.3,
    height= 11.7)
par(las= 2,
    mar= c(21,15,25,19),
    lwd= 0.5,
    cex.axis= 0.6)
plot(sel,
     padj_cutoff = 0.05,
     cex.balloons= 0.25)
title("GOs enrichment per cluster +/- PcG")
dev.off()



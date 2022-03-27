setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[project=="RNA_epiCancer"]

pdf("pdf/RNA/PCC_replicates.pdf")
meta[, {
    dds <- readRDS(dds_file)
    counts <- log2(counts(dds)+1)
    colnames(counts) <- gsub(".bam", "", colnames(counts))
    vl_heatmap(cor(counts),
               main= DESeq2_object)
    print("done")
}, .(dds_file, DESeq2_object)]
dev.off()

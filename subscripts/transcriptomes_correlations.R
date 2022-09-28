setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")

pdf("pdf/RNA_PCC_replicates.pdf")
par(las= 2,
    mar= c(7,7,5,6))
meta[, {
    dds <- readRDS(dds_file)
    counts <- log2(counts(dds)+1)
    colnames(counts) <- gsub(".bam", "", colnames(counts))
    vl_heatmap(cor(counts),
               main= DESeq2_object)
    print("done")
}, .(dds_file, DESeq2_object)]
dev.off()

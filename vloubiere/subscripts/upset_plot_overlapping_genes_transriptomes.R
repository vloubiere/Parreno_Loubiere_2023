setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[project=="RNA_epiCancer"]
meta <- meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), .(DESeq2_object, cdition)]
meta <- unique(meta[FC_file!="NA" & !grepl("vs_RNA_PHD11.txt$", FC_file)]) # Do not consider PHD11_T vs PHD11
meta <- meta[, fread(FC_file), (meta)]
meta[, cdition:= factor(cdition, 
                        levels= c("RNA_PH18",
                                  "RNA_PHD11",
                                  "RNA_PHD9",
                                  "RNA_PH29",
                                  "RNA_PHD11_T2",
                                  "RNA_PHD11_T3"))]
setorderv(meta, 
          c("DESeq2_object", 
            "cdition"))
meta <- meta[padj<0.05 & abs(log2FoldChange)>1]
meta[, class:= ifelse(log2FoldChange>0, "up", "down")]

pdf("pdf/RNA/upsetplot_overlapping_genes.pdf", width = 17)
layout(matrix(1:2, ncol= 2), 
       widths = c(1, 1.6))
meta[, {
    .l <<- split(FBgn, list(cdition, class))
    vl_upset_plot(.l, 
                  intersection_cutoff = 10)
}, DESeq2_object]
dev.off()
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
meta[, class:= ifelse(padj<0.05 & abs(log2FoldChange)>1, 
                      ifelse(log2FoldChange>0, "up", "down"), "unaffected")]
meta[is.na(padj), class:= "unaffected"]

pdf("pdf/Figures/alluvial_plot_timecourse_RNA.pdf", 
    width = 10, 
    height = 4)
layout(matrix(1:2, ncol= 2), 
       widths = c(1,1.3))
par(mar= c(1,1,3,1),
    cex= 0.7,
    las= 2)
meta[!grepl("vs_RNA_PHD11.txt$", FC_file), {
  .c <- dcast(.SD, FBgn~cdition, value.var = "class")[, -1]
  .c <- .c[apply(.c, 1, function(x) any(x!="unaffected"))]
  vl_alluvial_plot(.c,
                   class_levels = c("down", "unaffected", "up"),
                   col= c("cornflowerblue", "lightgrey", "tomato"))
  title(main= DESeq2_object, 
        line = 1.5,
        las= 1)
}, DESeq2_object]
dev.off()

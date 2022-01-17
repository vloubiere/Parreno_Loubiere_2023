setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

# FC
allograft <- data.table(file= c("db/FC_tables/epiCancer_ED_allograft_RNA_gDNA_RNA_PH18_vs_RNA_W18.txt",
                          "db/FC_tables/epiCancer_ED_allograft_RNA_gDNA_RNA_PHD11_vs_RNA_WKD.txt",
                          "db/FC_tables/epiCancer_ED_allograft_RNA_gDNA_RNA_PHD9_vs_RNA_WKD.txt",
                          "db/FC_tables/epiCancer_ED_allograft_RNA_gDNA_RNA_PH29_vs_RNA_W29.txt",
                          "db/FC_tables/epiCancer_ED_allograft_RNA_gDNA_RNA_PHD11_T2_vs_RNA_WKD.txt",
                          "db/FC_tables/epiCancer_ED_allograft_RNA_gDNA_RNA_PHD11_T3_vs_RNA_WKD.txt"),
                        cdition= c("PH18", "PHD11", "PHD9", "PH29", "PHD11_T2", "PHD11_T3"))
allograft <- allograft[, fread(file)[, .(FBgn, baseMean, log2FoldChange, padj)], cdition]
cutnrun_genotype <- data.table(file= c("db/FC_tables/epiCancer_ED_RNA_CUTNRUN_RNA_PH18_vs_RNA_W18.txt",
                                       "db/FC_tables/epiCancer_ED_RNA_CUTNRUN_RNA_PHD11_vs_RNA_WKD.txt",
                                       "db/FC_tables/epiCancer_ED_RNA_CUTNRUN_RNA_PHD9_vs_RNA_WKD.txt",
                                       "db/FC_tables/epiCancer_ED_RNA_CUTNRUN_RNA_PH29_vs_RNA_W29.txt"),
                               cdition= c("PH18", "PHD11", "PHD9", "PH29"))
cutnrun_genotype <- cutnrun_genotype[, fread(file)[, .(FBgn, baseMean, log2FoldChange, padj)], cdition]


# Object
allograft_dds <- readRDS("db/dds/epiCancer_ED_allograft_RNA_gDNA_dds.rds")
cutnrun_genotype_dds <- readRDS("db/dds/epiCancer_ED_RNA_CUTNRUN_dds.rds")
obj <- list(FC= list(allograft= allograft,
                     cutnrun_genotype= cutnrun_genotype),
            FPKM= list(allograft= as.data.table(DESeq2::fpkm(allograft_dds), keep.rownames= "FBgn"),
                       cutnrun= as.data.table(DESeq2::fpkm(cutnrun_genotype_dds), keep.rownames= "FBgn")),
            COUNTS= list(allograft= as.data.table(counts(allograft_dds), keep.rownames= "FBgn"),
                         cutnrun= as.data.table(counts(cutnrun_genotype_dds), keep.rownames= "FBgn")),
            VSD_COUNTS= list(allograft= as.data.table(SummarizedExperiment::assay(DESeq2::vst(allograft_dds, blind = FALSE)), keep.rownames= "FBgn"),
                             cutnrun= as.data.table(SummarizedExperiment::assay(DESeq2::vst(cutnrun_genotype_dds, blind = FALSE)), keep.rownames= "FBgn")))
# Rename columns
obj <- lapply(obj, function(x) lapply(x, function(y) setnames(y, gsub(".bam$", "", names(y)))))
saveRDS(obj, "Rdata/RNA_tables_object.rds")

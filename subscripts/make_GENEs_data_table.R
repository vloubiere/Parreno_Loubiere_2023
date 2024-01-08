setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import RNA-Seq FC ----
meta <- fread("Rdata/processed_metadata_RNA.txt")
meta <- meta[DESeq2_object=="epiCancer_noGFP" & FC_file!="NA"]
meta[, cdition:= factor(cdition, c("PH18", "PHD11", "PHD9", "PH29"))]
dat <- meta[, fread(FC_file), .(cdition, FC_file)]
dat <- dcast(dat, 
             FBgn~cdition, 
             value.var = list("log2FoldChange", "padj", "diff"))

# Add PHD11 vs PH29 comparison ----
PH29 <- fread("db/FC_tables/RNA/epiCancer_noGFP_PHD11_vs_PH29.txt")
setnames(PH29,
         c("FBgn", "log2FoldChange", "padj", "diff"),
         c("FBgn", "log2FoldChange_PHD11_vs_PH29", "padj_PHD11_vs_PH29", "diff_PHD11_vs_PH29"))
dat <- merge(dat, PH29, "FBgn")

# Add review RNA-Seq data ----
meta <- fread("Rdata/processed_metadata_RNA_review.txt")
meta <- meta[FC_file!="NA"]
meta <- unique(meta[, .(FC_file= unlist(tstrsplit(FC_file, ","))), cdition])
meta <- meta[grepl("_vs_W_|vs_GFP_W", FC_file)]
meta <- rbind(meta,
              data.table(cdition= c("ZFH1_PH_vs_GFP_PH", "STAT92E_PH_vs_GFP_PH"),
                         FC_file= c("db/FC_tables/RNA/epiCancer_rescue_ZFH1_PH_vs_GFP_PH.txt",
                                    "db/FC_tables/RNA/epiCancer_rescue_STAT92E_PH_vs_GFP_PH.txt")))
meta[, cdition:= factor(cdition, 
                        c("ZFH1_W", "STAT92E_W", "GFP_PH", "ZFH1_PH", "STAT92E_PH", "ZFH1_PH_vs_GFP_PH", "STAT92E_PH_vs_GFP_PH",
                          "PH_lL3",
                          "PH_eL3_0h_recov", "PH_eL3_24h_recov", "PH_eL3_48h_recov", "PH_eL3_96h_recov"))]
new <- meta[, fread(FC_file), .(cdition, FC_file)]
new <- dcast(new, 
             FBgn~cdition, 
             value.var = list("log2FoldChange", "padj", "diff"))
dat <- merge(dat, 
             new,
             by= "FBgn",
             all.x= T,
             all.y= T)

# Add RNA clusters ----
som <- readRDS("Rdata/clustering_RNA.rds")
dat[data.table(FBgn= rownames(som$data[[1]]), cl= som$unit.classif), cl:= i.cl, on= "FBgn"]
dat[, cl:= as.factor(cl)]
dat[, cl:= c("Transient-specific", "Irreversible", "Down 2", "Down 3", "Reversible", "Down 1")[as.numeric(cl)]]
dat[diff_PH29=="unaffected" & diff_PHD11=="unaffected" & diff_PHD9=="unaffected", cl:= "Unaffected"]
dat[, cl:= factor(cl, c("Unaffected", "Reversible", "Irreversible", "Transient-specific", "Down 1", "Down 2", "Down 3"))]

# Add FPKMs ----
dds <- readRDS("db/dds/RNA/epiCancer_noGFP_dds.rds")
fpkms <- as.data.table(DESeq2::fpkm(dds), keep.rownames = "FBgn")
fpkms <- melt(fpkms, id.vars = "FBgn")
fpkms[, variable:= tstrsplit(variable, "_", keep=1)]
fpkms <- fpkms[, .(FPKM= mean(value)), .(variable, FBgn)]
fpkms <- dcast(fpkms, FBgn~variable, value.var = "FPKM")
setnames(fpkms, names(fpkms)[-1], paste0("FPKM_", names(fpkms)[-1]))
dat <- fpkms[dat, on= "FBgn"]

# Add gene coordinates and symbols ----
gtf <- rtracklayer::import("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
gtf <- gtf[type=="gene", .(gene_id, symbol= gene_symbol, seqnames, start, end, strand)]
dat <- merge(dat,
             gtf,
             by.x= "FBgn",
             by.y= "gene_id")

# Add PcG binding ----
ov <- vl_closestBed(dat, "db/bed/merged_K27_domains/PH18.bed")
ov <- ov[dist==0]
ov[, max_start:= apply(.SD, 1, max), .SDcols= c("start", "start.b")]
ov[, min_end:= apply(.SD, 1, max), .SDcols= c("end", "end.b")]
ov <- ov[(min_end-max_start+1)/(end-start+1)>0.5] # 50% of the gene body covered with the mark
dat[, PcG_bound:= FBgn %in% ov$FBgn]

# Add RECOVERY ----
dat[, class:= fcase(PcG_bound & cl=="Irreversible", "Irreversible",
                    PcG_bound & cl=="Reversible", "Reversible",
                    PcG_bound & cl=="Unaffected", "Unaffected",
                    default= NA)]
dat[, class:= factor(class, c("Unaffected", "Irreversible", "Reversible"))]

# Save ----
setcolorder(dat,
            c("FBgn", "symbol", 
              "cl", "class", "PcG_bound",
              "seqnames", "start", "end", "strand", 
              grep("^diff", names(dat), value = T),
              grep("^log2FoldChange", names(dat), value = T),
              grep("^padj", names(dat), value = T)))
saveRDS(dat,
        "Rdata/final_gene_features_table.rds")
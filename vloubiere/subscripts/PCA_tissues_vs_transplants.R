setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import epiCancer FPKMs
dat <- lapply(c("db/dds/RNA/epiCancer_ED_allograft_RNA_gDNA_dds.rds",
                "db/dds/RNA/epiCancer_ED_RNA_CUTNRUN_dds.rds"), function(x) 
                {
                  dat <- DESeq2::fpkm(readRDS(x))
                  dat <- melt(as.data.table(dat, keep.rownames = T), id.vars = "rn")
                  dat[, variable:= gsub("_1.bam|_2.bam|_3.bam", "", variable)]
                  dat <- dat[, .(FPKM= mean(value)), .(rn, variable)]
                })
names(dat) <- c("GFP+_system", "GFP-_system")
dat <- rbindlist(dat, idcol = T)
dat[, variable:= gsub("^TPH", "PH", variable)]
dat[, variable:= gsub("^TW", "W", variable)]
dat[variable=="PHJ9", variable:= "PHD9"]
dat[variable=="PHJ11", variable:= "PHD11"]
dat[, col:= vl_palette_many_categ(12)[.GRP], variable]
dat <- dat[, .id:= paste0(variable, "_", .id)]
Cc <- unique(dat[, .(.id, col)])
dat <- dcast(dat, rn~.id, value.var = "FPKM")

# Import flybase tissues RNA-Seq
tissues <- fread("/mnt/d/_R_data/genomes/dm6/gene_rpkm_report_fb_2021_02.tsv")
tissues <- tissues[Parent_library_name=="modENCODE_mRNA-Seq_tissues"]
t.mat <- dcast(tissues, 
               `FBgn#`~RNASource_name, 
               value.var = "RPKM_value")
cols <- names(t.mat)[-1]
t.mat <- merge(t.mat, 
               dat,
               by.x= "FBgn#",
               by.y= "rn")
t.mat[is.na(t.mat)] <- 0
t.mat <- scale(as.matrix(t.mat, 1))
t.pca <- prcomp(t.mat)$rotation

plot(t.pca[,"PC1"],
     t.pca[,"PC2"],
     pch= NA)
text(t.pca[,"PC1"],
     t.pca[,"PC2"],
     gsub("mE_mRNA_|_system", "", rownames(t.pca)),
     cex= 1,
     col= ifelse(grepl("_system$", rownames(t.pca)), "red", "black"))



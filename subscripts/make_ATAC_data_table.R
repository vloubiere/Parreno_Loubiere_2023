setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)

# Import all tumor peaks ----
dat <- fread("db/saf/ATAC_peaks_tumor.saf", col.names = c("ID", "seqnames", "start", "end"), sel= 1:4)

# Add peaks overlaps ----
# ATAC-Seq
peaks <- data.table(file= list.files("db/peaks/ATAC/", "^ATAC_PH.*conf", full.names = T))
peaks[, cdition:= paste0(gsub("^ATAC_|_conf_peaks.narrowPeak$", "", basename(file)), "_ATAC")]
peaks <- peaks[cdition %in% c("PH18_ATAC", "PH29_ATAC", "PHD11_ATAC")]
peaks <- peaks[, .(ID= dat$ID, vl_covBed(dat, file)>0), cdition]
peaks <- dcast(peaks, ID~cdition, value.var = "V2")
dat <- merge(dat,
             peaks,
             by= "ID",
             all.x= T)

# Add clusters ----
som <- readRDS("Rdata/clustering_ATAC.rds")
dat[rownames(som$data[[1]]), som:= som$unit.classif, on= "ID"]
dat[is.na(som), som:= 0]
.n <- c("Unaffected", "Irreversible", "Reversible", "Decreased")
dat[, cl:= .n[som+1]]
dat[, cl:= factor(cl, .n)]
dat$som <- NULL

# Add ATAC-Seq FC ----
# Tumors 
FC <- data.table(FC_file= list.files("db/FC_tables/ATAC/", "^ATAC.*_vs_PH18|^ATAC.*_vs_PH29", full.names = T))
FC[, cdition:= gsub("^ATAC_|_vs_PH18|.txt$", "", basename(FC_file))]
FC <- FC[, fread(FC_file), cdition]
FC[, diff:= fcase(log2FoldChange>1 & padj<0.001, "up",
                  log2FoldChange<(-1) & padj<0.001, "down",
                  default = "unaffected")]
FC <- dcast(FC,
            ID+baseMean~cdition,
            value.var = list("log2FoldChange", "padj", "diff"))
dat <- merge(dat,
             FC,
             by= "ID",
             all.x= T)

# Save ----
saveRDS(dat,
        "Rdata/final_ATAC_table.rds")
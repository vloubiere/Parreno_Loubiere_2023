setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(glmnet)

# Import peaks ----
dat <- readRDS("Rdata/final_ATAC_rescue_table.rds")
peaks <- vl_resizeBed(dat, "center", 0, 0)
peaks <- peaks[diff_zfh1_ph=="unaffected"]
peaks <- peaks[diff_zfh1_ph_vs_gfp_ph=="up" & ZFH1_counts>0]

# Import genes ----
genes <- readRDS("Rdata/final_gene_features_table.rds")
TSS <- vl_resizeBed(genes, "start", 0, 0)
cl <- vl_closestBed(peaks, TSS)
cl <- cl[abs(dist)<25e3]

# GO enrichment ----
# enr <- vl_GO_enrich(cl[class=="Increased upon zfh1-KD", FBgn.b],
enr <- vl_GO_enrich(cl$FBgn.b,
                    geneUniverse.IDs = TSS$FBgn,
                    species = "Dm",
                    select = "BP")

# Plot ----
pdf("pdf/review_GO_enrichment_rescued_peaks.pdf", 4, 3)
vl_par(mai= c(.9,2,.9,1.5))
plot(enr,
     top.enrich= 8,
     padj.cutoff= 0.05,
     min.counts= 10,
     order= "log2OR",
     col= c("cornflowerblue", "tomato"))
dev.off()
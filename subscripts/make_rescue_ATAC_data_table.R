setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import FC
dat <- list("gfp_ph"= "db/FC_tables/ATAC/ATAC_GFP_PH_vs_GFP_W.txt",
            "stat_ph"= "db/FC_tables/ATAC/ATAC_Stat_PH_vs_GFP_W.txt",
            "zfh1_ph"= "db/FC_tables/ATAC/ATAC_ZFH1_PH_vs_GFP_W.txt",
            "stat_ph_vs_gfp_ph"= "db/FC_tables/ATAC/ATAC_Stat_PH_vs_GFP_PH.txt",
            "zfh1_ph_vs_gfp_ph"= "db/FC_tables/ATAC/ATAC_ZFH1_PH_vs_GFP_PH.txt")
dat <- lapply(dat, fread)
dat <- rbindlist(dat, idcol = "variable")
# Define diff peaks
dat[, diff:= fcase(log2FoldChange>1 & padj<0.001, "up",
                   log2FoldChange<(-1) & padj<0.001, "down",
                   default= "unaffected")]
# Cast
dat <- dcast(dat,
             ID+baseMean~variable,
             value.var = c("log2FoldChange", "padj", "diff"))
dat[, c("seqnames", "start", "end"):= vl_toDTranges(ID)]
# Add zfh1 and stat peaks
res <- vl_resizeBed(dat,
                    "center",
                    250,
                    250,
                    genome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
counts <- vl_motif_counts(res,
                          sel= c("cisbp__M4483", "jaspar__MA0532.1"),
                          genome = BSgenome.Dmelanogaster.UCSC.dm6::BSgenome.Dmelanogaster.UCSC.dm6)
setnames(counts, c("ZFH1_counts", "STAT92E_counts"))
dat <- cbind(dat, counts)

saveRDS(dat, "Rdata/final_ATAC_rescue_table.rds")
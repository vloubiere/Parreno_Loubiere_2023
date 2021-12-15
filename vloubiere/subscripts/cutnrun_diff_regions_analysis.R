setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(rtracklayer)
require(GenomicRanges)

# Cutnrun FC
FC <- fread("Rdata/K27_cutnrun_FC.txt")
FC$H3K27Ac_PH18 <- log2(vl_bw_coverage(FC, "db/bw/cutnrun_merge/H3K27Ac_PH18_merge.bw"))
FC$H3K27me3_PH18 <- log2(vl_bw_coverage(FC, "db/bw/cutnrun_merge/H3K27me3_PH18_merge.bw"))
FC$PH <- log2(vl_bw_coverage(FC, "db/bw/SA_2020/PH_ED_merge.bw"))
FC$SUZ12 <- log2(vl_bw_coverage(FC, "db/bw/SA_2020/SUZ12_ED_merge.bw"))

# Genes FC
RNA <- readRDS("Rdata/final_FC_table.rds")
genes_coor <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes_coor) <- "UCSC"
genes_coor <- as.data.table(genes_coor)[type%in% c("mRNA", "gene")]
RNA <- genes_coor[RNA, unique(data.table(seqnames, start, end, strand, cdition= i.cdition, FBgn, transcript_id, symbol, log2FoldChange)), on= "gene_id==FBgn", allow.cartesian= TRUE]
RNA[strand=="+", c("start", "end"):= .(start-5000, start+2000)]
RNA[strand=="-", c("start", "end"):= .(end-2000, end+5000)]
RNA <- dcast(RNA[grep("^PH.*TS$", cdition)], 
             seqnames+start+end+FBgn+transcript_id~cdition, 
             value.var= "log2FoldChange")

# Add RNA FC to ChIP
res <- cbind(FC, RNA[FC, .(FBgn= .(FBgn[!duplicated(FBgn)]), 
                           PH18_TS= mean(PH18_TS[!duplicated(FBgn)], na.rm=T), 
                           PH29_TS= mean(PH29_TS[!duplicated(FBgn)], na.rm=T), 
                           PHD9_TS= mean(PHD9_TS[!duplicated(FBgn)], na.rm=T), 
                           PHD11_TS= mean(PHD11_TS[!duplicated(FBgn)], na.rm=T)), 
                     .EACHI, 
                     on= c("seqnames", "start<end", "end>end")][, .(FBgn, PH18_TS, PH29_TS, PHD9_TS, PHD11_TS)])
res[, coor:= paste0(seqnames, ":", start, "-", end)]
K27me3 <- res[region_type=="K27me3" & (abs(H3K27me3_PH29_log2FC)>1 | abs(H3K27me3_PHD9_log2FC)>1 | abs(H3K27me3_PHD11_log2FC)>1)]
K27Ac <- res[region_type=="K27Ac" & (abs(H3K27Ac_PH29_log2FC)>1 | abs(H3K27Ac_PHD9_log2FC)>1 | abs(H3K27Ac_PHD11_log2FC)>1)]

pdf("pdf/cutnrun/heatmap_diff_regions_cutnrun_K27.pdf", 
    height= 5)
par(mfrow= c(1,4))
#--------------------------#
# K27me3 diff regions
#--------------------------#
# K27me FC
K27me3_cl <- vl_heatmap(K27me3[, .(coor, H3K27me3_PH29_log2FC, H3K27me3_PHD11_log2FC, H3K27me3_PHD9_log2FC)], 
                 kmeans= 4,
                 breaks= c(-2.5,-0.5,0.5,2.5),
                 main= "K27me3 diff.",
                 legend_title= "K27me3 log2FC", 
                 show_rownames= F,
                 col= c("cornflowerblue", "white", "white", "tomato"))
plot(K27me3_cl,
     show_rownames= F,
     main= "WT enrich.",
     add= as.matrix(K27me3[, .(coor, 
                               H3K27Ac_PH18, 
                               H3K27me3_PH18, 
                               SUZ12_WT= SUZ12, 
                               PH_WT= PH)], 1),
     legend_title= "log2 enrich.",
     auto_margins= F)
# K27Ac FC
plot(K27me3_cl,
     main= "K27Ac FC",
     as.matrix(K27me3[, .(coor, H3K27Ac_PH29_log2FC, H3K27Ac_PHD11_log2FC, H3K27Ac_PHD9_log2FC)], 1),
     show_rownames= F,
     breaks= c(-2.5,-0.5,0.5,2.5), 
     auto_margins= F,
     legend_title= "K27Ac log2FC", 
     show_col_dendrogram= F,
     col= c("cornflowerblue", "white", "white", "tomato"))
mat_FC <- as.matrix(K27me3[, .(coor, PH18_TS, PH29_TS, PHD9_TS, PHD11_TS)], 1)
mat_FC[is.nan(mat_FC)] <- NA
plot(K27me3_cl,
     main= "Transcript. FC",
     mat_FC,
     show_rownames= F, 
     breaks= c(-1,-0.25,0.25,1), 
     auto_margins= F,
     legend_title= "RNA log2FC", 
     show_col_dendrogram= F,
     col= c("cornflowerblue", "white", "white", "tomato"))

#--------------------------#
# K27Ac diff regions
#--------------------------#
# K27Ac FC
K27Ac_cl <- vl_heatmap(K27Ac[, .(coor, H3K27Ac_PH29_log2FC, H3K27Ac_PHD11_log2FC, H3K27Ac_PHD9_log2FC)], 
                       kmeans= 4,
                       breaks= c(-2.5,-0.5,0.5,2.5),
                       main= "K27Ac diff.",
                       legend_title= "K27Ac log2FC", 
                       show_rownames= F,
                       col= c("cornflowerblue", "white", "white", "tomato"))
plot(K27Ac_cl,
     show_rownames= F,
     main= "WT enrich.",
     add= as.matrix(K27Ac[, .(coor, 
                               H3K27Ac_PH18= H3K27Ac_PH18, 
                               H3K27me3_PH18= H3K27me3_PH18, 
                               SUZ12_WT= SUZ12, 
                               PH_WT= PH)], 1),
     legend_title= "log2 enrich.",
     auto_margins= F)
# K27me3 FC
plot(K27Ac_cl,
     main= "K27me3 FC",
     as.matrix(K27Ac[, .(coor, H3K27me3_PH29_log2FC, H3K27me3_PHD11_log2FC, H3K27me3_PHD9_log2FC)], 1),
     show_rownames= F,
     breaks= c(-2.5,-0.5,0.5,2.5), 
     auto_margins= F,
     legend_title= "K27me3 log2FC", 
     show_col_dendrogram= F,
     col= c("cornflowerblue", "white", "white", "tomato"))
# Transcript FC
mat_FC <- as.matrix(K27Ac[, .(coor, PH18_TS, PH29_TS, PHD9_TS, PHD11_TS)], 1)
mat_FC[is.nan(mat_FC)] <- NA
plot(K27Ac_cl,
     main= "Transcript. FC",
     mat_FC,
     show_rownames= F, 
     breaks= c(-1,-0.25,0.25,1), 
     auto_margins= F,
     legend_title= "RNA log2FC", 
     show_col_dendrogram= F,
     col= c("cornflowerblue", "white", "white", "tomato"))
dev.off()

saveRDS(unique(K27Ac_cl$result_DT),
        "Rdata/K27Ac_cutnrun_clustering.rds")
saveRDS(unique(K27Ac_cl$result_DT[, .(rn, rn_cl)]),
        "Rdata/K27Ac_cutnrun_clustering_simplif.rds")
saveRDS(unique(K27me3_cl$result_DT),
        "Rdata/K27me3_cutnrun_clustering.rds")
saveRDS(unique(K27Ac_cl$result_DT[, .(rn, rn_cl)]),
        "Rdata/K27Ac_cutnrun_clustering_simplif.rds")

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)

# Import clusters
K27me3 <- readRDS("Rdata/K27me3_cutnrun_clustering.rds")
K27me3 <- cbind(as.data.table(GRanges(K27me3$rn)), K27me3[,"rn_cl"])

K27Ac <- readRDS("Rdata/K27Ac_cutnrun_clustering.rds")
K27Ac <- cbind(as.data.table(GRanges(K27Ac$rn)), K27Ac[,"rn_cl"])

# Import genes coor
genes_coor <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes_coor) <- "UCSC"
genes_coor <- as.data.table(genes_coor)[type=="mRNA"]
genes_coor[strand=="+", c("start", "end"):= .(start-5000, start+2000)]
genes_coor[strand=="-", c("start", "end"):= .(end-2000, end+5000)]

# Overlap
res_K27me3 <- genes_coor[K27me3, 
                         .(i.rn_cl, FBgn= gene_id, gene_symbol), 
                         .EACHI, 
                         on= c("seqnames", "end>start", "start<start")]
res_K27me3 <- na.omit(unique(res_K27me3[, .(FBgn, cl= i.rn_cl)]))
res_K27Ac <- genes_coor[K27Ac, 
                        .(i.rn_cl, FBgn= gene_id, gene_symbol), 
                        .EACHI, 
                        on= c("seqnames", "end>start", "start<start")]
res_K27Ac <- na.omit(unique(res_K27Ac[, .(FBgn, cl= i.rn_cl)]))


pdf("pdf/cutnrun/GO_K27_cluster_genes.pdf", 
    width = 10, height = 10)
vl_GO_clusters(split(res_K27Ac$FBgn, res_K27Ac$cl), 
               padj_cutoff = 0.00001, 
               N_top = 10,
               cex.balloons = 1.5, 
               main = "K27Ac")
vl_GO_clusters(FBgn_list = split(res_K27me3$FBgn, res_K27me3$cl), 
               padj_cutoff = 0.01,
               N_top = 10,
               cex.balloons = 1.5, 
               auto_margins = F, 
               main = "K27me3")
dev.off()

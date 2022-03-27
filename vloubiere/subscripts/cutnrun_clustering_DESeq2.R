setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(kohonen)

# Import
FBGN <- rtracklayer::import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
GenomeInfoDb::seqlevelsStyle(FBGN) <- "UCSC"
FBGN <- as.data.table(FBGN)
promoters <- FBGN[type=="gene", .(seqnames, 
                                  start= ifelse(strand=="+", start, end), 
                                  strand), .(gene_id, gene_symbol)]
promoters[, end:= start]
meta <- fread("Rdata/processed_metadata_CUTNRUN.txt")

# Import mut FC and select diff regions
mut <- meta[!is.na(FC_file_K27_segmentation), {
  fread(.BY[[1]], colClasses = c("character", rep("numeric", 5)))
}, .(FC_file_K27_segmentation, ChIP, cdition)]
mut[is.na(padj) | padj>0.05, log2FoldChange:= NA]
# mut <- mut[mut[, any(!is.na(log2FoldChange)), coor][(V1), coor], on= "coor"]

# Import WT levels (ChIP vs IgG PH18)
WT <- vl_importBed("db/narrowpeaks/cutnrun/K27_segmentation_mergedBed.bed")
WT[, coor:= paste0(seqnames, ":", start, "-", end)]

# Format data
dat <- dcast(mut, coor~ChIP+cdition, value.var = "log2FoldChange")
dat <- merge(WT[, .(coor, K27me3_log2_enr_PH18, K27Ac_log2_enr_PH18)], dat, by= "coor")

# Clustering
layers_NA <- list(WT= c("coor", grep("log2_enr", names(dat), value = T)),
                  constant= c("coor", grep("^H3.*_PH29$", names(dat), value = T)),
                  transient= c("coor", grep("^H3.*_PHD.*$", names(dat), value = T)))
# layers <- lapply(layers, function(x) scale(as.matrix(dat[, ..x], 1)))
layers_NA <- lapply(layers_NA, function(x) as.matrix(dat[, ..x], 1))
layers_NA[[1]] <- scale(layers_NA[[1]])
layers <- lapply(layers_NA, function(x) {x[is.na(x)] <- 0; return(x)})
grid <- somgrid(2, 
                3, 
                "hexagonal", 
                toroidal= T)
init <- lapply(layers, function(x)
{
  set.seed(6)
  x[sample(nrow(x), grid$xdim*grid$ydim),]
})
som <- supersom(data = layers, 
                grid= grid,
                init = init,
                user.weights= c(1,1,3), 
                maxNA.fraction = 1)
som$ordered_matrices <- do.call(cbind, layers_NA)[order(som$unit.classif),]

# Cluster table
cl <- as.data.table(GRanges(dat$coor))
cl[, cl:= som$unit.classif]
setorderv(cl, "cl")

# Closest genes
genes <- vl_closestBed(cl,
                       promoters)
genes <- genes[, .SD[abs(dist)==min(abs(dist)),
                     .(FBgn= gene_id.b, 
                       symbol= gene_symbol.b, 
                       dist)], .(seqnames, start= start.a, end= end.a, cl= cl.a)]

# Final object
obj <- list(regions= cl,
            genes= genes,
            som= som,
            plot_cluster= function(som)
            {
              # Plot
              cuts <- nrow(som$ordered_matrices)-cumsum(table(som$unit.classif))+0.5
              Cc <- c("cornflowerblue", "white", "tomato")
              
              layout(matrix(1:2, ncol= 2),
                     widths = c(1,1.45))
              par(mar= c(11.5, 4.5, 1.5, 7))
              vl_heatmap(som$ordered_matrices[, 1:2],
                         cluster_rows= F,
                         cluster_cols = F, 
                         breaks = c(-2, 0, 2), 
                         col = Cc,
                         legend_title = "PH18 enr. (log2)", 
                         auto_margins = F, 
                         show_rownames = F)
              cl_vec <- length(som$unit.classif)-(cumsum(table(som$unit.classif))-table(som$unit.classif)/2)
              axis(2, 
                   at= cl_vec, 
                   labels = names(cl_vec), 
                   las= 1, lwd= NA)
              abline(h= cuts)
              par(mar= c(11.5, 0, 1.5, 8))
              vl_heatmap(som$ordered_matrices[, 3:8],
                         cluster_rows= F,
                         cluster_cols = F, 
                         show_rownames = F, 
                         auto_margins = F, 
                         breaks = c(-1.5, 0, 1.5),
                         col= Cc,
                         legend_title = "FC over PH18 (log2)")
              abline(h= cuts)
              abline(v= 2.5)
            })

saveRDS(obj,
        "Rdata/clustering_cutnrun.rds")

clusters <- split(obj$regions, obj$regions$cl)
lapply(seq(clusters), function(i) vl_exportBed(clusters[[i]][, 1:3], paste0("db/peaks/K27_cutnrun_changes/cluster_", i, "_DESeq2.bed")))

pdf("pdf/cutnrun/clustering_diff_regions.pdf")
obj$plot_cluster(obj$som)
dev.off()





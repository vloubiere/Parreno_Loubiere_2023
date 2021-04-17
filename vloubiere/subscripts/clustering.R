# Import conditions of interest
dat <- data.table(files= sapply(c("EzJ9_vs_WKD", 
                                  "EzJ11_vs_WKD", 
                                  "Ez29_vs_W29", 
                                  "PHJ9_vs_WKD", 
                                  "PHJ11_vs_WKD", 
                                  "PH29_vs_W29"), function(x) 
                                    list.files("db/FC_tables/", x, full.names = T)))
dat[, cdition:= tstrsplit(files, "epiCancer_|_FC", keep= 2)]
dat <- dat[, fread(files, 
                   select = c(1,2,3,6), 
                   col.names = c("FBgn", "baseMean", "log2FoldChange", "padj")),  (dat)]

# Select genes differentially expressed in at least one of the 6 conditions
sel <- dat[, any(padj<0.001 & abs(log2FoldChange)>log2(1.5)), FBgn][(V1), FBgn]
diff <- dat[FBgn %in% sel]

# Remove genes deferentially expressed at 18
sel <- rbind(fread("/_R_data/projects/epigenetic_cancer/db/FC_tables/RNA_epiCancer_PH18_vs_W18_FC.txt"),
             fread("/_R_data/projects/epigenetic_cancer/db/FC_tables/RNA_epiCancer_Ez18_vs_W18_FC.txt"))
sel <- sel[padj<=0.05 & abs(log2FoldChange)>log2(1.5), V1]
clean <- diff[!(FBgn %in% sel)]

# Clip extremes
vl_boxplot(clean, log2FoldChange~cdition, outline = T)
abline(h= c(-8, 8))
clean[, clipped_log2FoldChange:= log2FoldChange]
clean[log2FoldChange> 8, clipped_log2FoldChange:= 8]
clean[log2FoldChange< -8, clipped_log2FoldChange:= -8]

# Clustering
grid <- somgrid(4, 4, "hexagonal", toroidal= T)
train <- as.matrix(dcast(clean, FBgn~cdition, value.var = "clipped_log2FoldChange"), 1)
set.seed(5)
som <- supersom(train, grid= grid, maxNA.fraction = 0.1, user.weights = c(1, 10))

# Add cluster columns to the data
som_cl <- as.data.table(som$data, keep.rownames = "FBgn")
som_cl[, cl:= som$unit.classif]
clean[som_cl, cl:= i.cl, on= "FBgn"]
clean[, cdition:= factor(cdition, 
                         levels = c("EzJ9_vs_WKD", 
                                    "EzJ11_vs_WKD", 
                                    "Ez29_vs_W29", 
                                    "PHJ9_vs_WKD", 
                                    "PHJ11_vs_WKD", 
                                    "PH29_vs_W29"))]
setorderv(clean, c("cl", "FBgn"))

# Add genes symbols
symbols <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
clean[symbols, symbol:= i.gene_symbol, on= "FBgn==gene_id"]

# Save tables and som
saveRDS(som, "Rdata/som_clustering_transcriptomes.rds")
saveRDS(clean, "Rdata/final_clustering_transcriptomes.rds")

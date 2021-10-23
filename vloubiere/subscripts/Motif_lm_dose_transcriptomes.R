setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
 
FC <- readRDS("Rdata/final_FC_table.rds")
feat <- readRDS("Rdata/ED_REs_features_final.rds")
if(!exists("mot"))
  mot <- readRDS("Rdata/ED_REs_motif_counts_final.rds")
dat <- data.table(file= list.files("db/motif_enrichment_tables/",
                                   "dose_clusters", 
                                   full.names = T))
dat <- dat[, readRDS(file), file]
dat <- mot[motif %in% dat[padj<0.001, motif]]
# Select only open regions?
# dat[feat, open:= i.open, on= "RE_ID"]
# dat <- dat[, .(motif_counts= max(c(motif_counts[open], 0))), .(FBgn, motif, motif_name)]
dat <- dat[, .(motif_counts= max(motif_counts)), .(FBgn, motif, motif_name)]
dat[FC[cdition=="PH29_dose"], c("log2FoldChange", "padj"):= .(i.log2FoldChange, i.padj), on= "FBgn"]
dat <- na.omit(dat)
dat[, motif_counts:= ifelse(motif_counts>quantile(motif_counts, 0.98), 
                            quantile(motif_counts, 0.98), 
                            motif_counts), motif]
# test <- dat[, summary(lm(log2FoldChange~motif_counts, .SD))$r.squared, .(motif, motif_name)]
test <- dat[, lm(log2FoldChange~motif_counts, .SD)$coefficients[2], .(motif, motif_name)]
test <- test[, .SD[order(abs(V1), decreasing = T)][1], motif_name]

dir.create("pdf/modelling/", 
           showWarnings = F)

pdf("pdf/modelling/boxplots_lm_top_motifs.pdf", 
    width = 16, 
    height = 5)
par(mfrow= c(2, 10))
setorderv(test, "V1", order = -1)
test[1:10, {
  if(V1>0)
    boxplot(log2FoldChange~motif_counts, 
            dat[.BY, , on= "motif"], 
            outline= F, 
            main= motif_name)
  else
    plot.new()
  print("")
}, .(motif, motif_name)]
setorderv(test, "V1", order = 1)
test[1:10, {
  boxplot(log2FoldChange~motif_counts, 
          dat[.BY, , on= "motif"], 
          outline= F, 
          main= motif_name)
  print("")
}, .(motif, motif_name)]
dev.off()


dat <- data.table(file= list.files("db/counts/RNA_epiCancer/", "^PH29|^PHJ11|^PHD11", full.names = T))
dat <- dat[, {
  x <- as.data.table(readRDS(file)$counts, keep.rownames= "FBgn")
  colnames(x)[2] <- "counts"
  x
}, file]
dat[, cdition:= gsub("_counts.rds", "", basename(file)), file]
dat <- dcast(dat, FBgn~cdition, value.var = "counts")
dat <- dat[, .(FBgn, sapply(.SD, function(x) {
  x <- x+1
  return(x/sum(x)*1e6)
})), .SDcols= patterns("^PH")]


pdf("pdf/PCC_transplant_vs_origin.pdf")
par(mar= c(9,9,4,5))
vl_heatmap(cor(dat[,!"FBgn"]), display_numbers = T, col= c("blue", "yellow"))
dev.off()

# Compute raw log2FC
diff <- list(PHD11_T5_1_2= data.table(FBgn= dat$FBgn, 
                           log2FC= rowMeans(dat[, .(log2(PHD11_T5_1_2)-log2(PHJ11_1), log2(PHD11_T5_1_2)-log2(PHJ11_2), log2(PHD11_T5_1_2)-log2(PHJ11_3))])),
             PHD11_T8_1_2= data.table(FBgn= dat$FBgn, 
                                      log2FC= rowMeans(dat[, .(log2(PHD11_T8_1_2_1)-log2(PHJ11_1), log2(PHD11_T8_1_2_1)-log2(PHJ11_2), log2(PHD11_T8_1_2_1)-log2(PHJ11_3))])),
             PH29_T5_1_1= data.table(FBgn= dat$FBgn, 
                                     log2FC= rowMeans(dat[, .(log2(PH29_T5_1_1)-log2(PH29_1), log2(PH29_T5_1_1)-log2(PH29_2), log2(PH29_T5_1_1)-log2(PH29_3))])),
             PH29_T5_2_3= data.table(FBgn= dat$FBgn, 
                                     log2FC= rowMeans(dat[, .(log2(PH29_T5_2_3)-log2(PH29_1), log2(PH29_T5_2_3)-log2(PH29_2), log2(PH29_T5_2_3)-log2(PH29_3))])),
             PH29_T8_2_3= data.table(FBgn= dat$FBgn, 
                                     log2FC= rowMeans(dat[, .(log2(PH29_T8_2_3_1)-log2(PH29_1), log2(PH29_T8_2_3_1)-log2(PH29_2), log2(PH29_T8_2_3_1)-log2(PH29_3))])))
diff <- rbindlist(diff, idcol = "cdition")
diff <- dcast(diff, FBgn~cdition, value.var = "log2FC")
diff <- diff[apply(diff[, !"FBgn"], 1, function(x) any(abs(x)>2))]
symbols <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
diff[symbols, symbol:= i.gene_symbol, on= "FBgn==gene_id"]

# # Clustering
# mat <- as.matrix(diff[, -c("FBgn", "symbol")], paste0(diff$FBgn, "__", diff$symbol))
# pdf("pdf/clustering_transplantation_transcriptomes.pdf", height= 12)
# par(mar= c(9,12,2,5))
# Tcl <- vl_heatmap(mat, cutree_rows = 8)
# dev.off()
# saveRDS(Tcl, "Rdata/clustering_transplantations.rds")

cl <- readRDS("Rdata/clustering_transplantations.rds")
boxplot(value~rcl, cl)





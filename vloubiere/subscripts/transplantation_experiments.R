# DESeq2 transplantations
counts <- data.table(file= list.files("db/counts/RNA_epiCancer/", "W29|PH29_T|WKD|PHD11_T", full.names = T))
counts[grepl("PH29_T", basename(file)), cdition:= "PH29T"]
counts[grepl("W29", basename(file)), cdition:= "W29"]
counts[grepl("PHD11_T", basename(file)), cdition:= "PHD11T"]
counts[grepl("WKD", basename(file)), cdition:= "WKD"]
counts[, rep:= seq(.N), cdition]
counts <- counts[, as.data.table(readRDS(file)$counts, keep.rownames= "FBgn"), .(file, cdition, rep)]
colnames(counts)[ncol(counts)] <- "counts"
DF <- data.frame(dcast(counts, FBgn~cdition+rep, value.var = "counts"), row.names = 1)
sampleTable <- data.frame(unique(counts[, .(cdition, rep)]))
rownames(sampleTable) <- paste0(sampleTable$cdition, "_", sampleTable$rep)
sampleTable <- sampleTable[match(colnames(DF),rownames(sampleTable)),]
DF <- DF[rowSums(DF)>10,]
dds <- DESeqDataSetFromMatrix(countData= DF, 
                              colData= sampleTable, 
                              design= ~rep+cdition)
dds <- DESeq(dds)
saveRDS(dds, "db/dds/RNA_transplantations_dds.rds")


.c <- lfcShrink(dds,
                type= "ashr",
                contrast= c("cdition", "PHD11T", "WKD"))
.c <- as.data.table(as.data.frame(.c), keep.rownames = "V1")
fwrite(.c, 
       file = "db/FC_tables/RNA_epiCancer_PHD11T_vs_WKD_FC.txt", 
       col.names = T, 
       row.names = F, 
       sep= "\t", 
       quote= F)
.c <- lfcShrink(dds,
                type= "ashr",
                contrast= c("cdition", "PH29T", "W29"))
.c <- as.data.table(as.data.frame(.c), keep.rownames = "V1")
fwrite(.c, 
       file = "db/FC_tables/RNA_epiCancer_PH29T_vs_W29_FC.txt", 
       col.names = T, 
       row.names = F, 
       sep= "\t", 
       quote= F)

# # RAW logFC
# dat <- data.table(file= list.files("db/counts/RNA_epiCancer/", "^PH29|^PHJ11|^PHD11", full.names = T))
# dat <- dat[, {
#   x <- as.data.table(readRDS(file)$counts, keep.rownames= "FBgn")
#   colnames(x)[2] <- "counts"
#   x
# }, file]
# dat[, cdition:= gsub("_counts.rds", "", basename(file)), file]
# dat <- dcast(dat, FBgn~cdition, value.var = "counts")
# dat <- dat[, .(FBgn, sapply(.SD, function(x) {
#   x <- x+1
#   return(x/sum(x)*1e6)
# })), .SDcols= patterns("^PH")]
# 
# 
# pdf("pdf/PCC_transplant_vs_origin.pdf")
# par(mar= c(9,9,4,5))
# vl_heatmap(cor(dat[,!"FBgn"]), display_numbers = T, col= c("blue", "yellow"))
# dev.off()
# 
# # Compute raw log2FC
# diff <- list(PHD11_T5_1_2= data.table(FBgn= dat$FBgn, 
#                            log2FC= rowMeans(dat[, .(log2(PHD11_T5_1_2)-log2(PHJ11_1), log2(PHD11_T5_1_2)-log2(PHJ11_2), log2(PHD11_T5_1_2)-log2(PHJ11_3))])),
#              PHD11_T8_1_2= data.table(FBgn= dat$FBgn, 
#                                       log2FC= rowMeans(dat[, .(log2(PHD11_T8_1_2_1)-log2(PHJ11_1), log2(PHD11_T8_1_2_1)-log2(PHJ11_2), log2(PHD11_T8_1_2_1)-log2(PHJ11_3))])),
#              PH29_T5_1_1= data.table(FBgn= dat$FBgn, 
#                                      log2FC= rowMeans(dat[, .(log2(PH29_T5_1_1)-log2(PH29_1), log2(PH29_T5_1_1)-log2(PH29_2), log2(PH29_T5_1_1)-log2(PH29_3))])),
#              PH29_T5_2_3= data.table(FBgn= dat$FBgn, 
#                                      log2FC= rowMeans(dat[, .(log2(PH29_T5_2_3)-log2(PH29_1), log2(PH29_T5_2_3)-log2(PH29_2), log2(PH29_T5_2_3)-log2(PH29_3))])),
#              PH29_T8_2_3= data.table(FBgn= dat$FBgn, 
#                                      log2FC= rowMeans(dat[, .(log2(PH29_T8_2_3_1)-log2(PH29_1), log2(PH29_T8_2_3_1)-log2(PH29_2), log2(PH29_T8_2_3_1)-log2(PH29_3))])))
# diff <- rbindlist(diff, idcol = "cdition")
# diff <- dcast(diff, FBgn~cdition, value.var = "log2FC")
# diff <- diff[apply(diff[, !"FBgn"], 1, function(x) any(abs(x)>2))]
# symbols <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
# diff[symbols, symbol:= i.gene_symbol, on= "FBgn==gene_id"]
# 
# # # Clustering
# # mat <- as.matrix(diff[, -c("FBgn", "symbol")], paste0(diff$FBgn, "__", diff$symbol))
# # pdf("pdf/clustering_transplantation_transcriptomes.pdf", height= 12)
# # par(mar= c(9,12,2,5))
# # Tcl <- vl_heatmap(mat, cutree_rows = 8)
# # dev.off()
# # saveRDS(Tcl, "Rdata/clustering_transplantations.rds")
# 
# cl <- readRDS("Rdata/clustering_transplantations.rds")
# boxplot(value~rcl, cl)




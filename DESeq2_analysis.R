setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(DESeq2)
require(pheatmap)
require(data.table)
require(gridExtra)
require(clusterProfiler)
require(org.Dm.eg.db)
require(rtracklayer)
require(colorRamps)

#--------------------------------------------------------------------#
# 1- Alignment statistics
#--------------------------------------------------------------------#
stat <- fread("alignment_stats/STATISTICS_TABLE.txt")
setcolorder(stat, c(1,2,4,3))
cols <- colnames(stat)[2:3]
stat[, (cols):= lapply(.SD, function(x) formatC(x, big.mark = ",")), .SDcols= cols]

pdf("pdf/alignment_statistics.pdf", width = 9)
grid.table(stat)
dev.off()

#--------------------------------------------------------------------#
# 2- DESeq2 analysis
#--------------------------------------------------------------------#
# Import and clean counts
dat <- data.table(file= list.files("counts", ".txt", full.names = T))
dat[, cdition:= gsub("_counts.txt", "", basename(file))]
dat <- dat[, fread(file), .(cdition, file)]
colnames(dat)[3:4] <- c("symbol", "counts")
dat <- dat[!grepl("^__", symbol)]

# Correlations
dmat <- dcast(dat, symbol~cdition, value.var= "counts")
mat <- as.matrix(dmat, 1)
pheatmap(cor(mat), filename = "pdf/PCC_counts.pdf")

# DESeq2 analysis
if(!file.exists("data/dds_object_1.0.rds"))
{
  DF <- as.data.frame(mat)
  DF <- DF[rowSums(DF)>5,]
  sampleTable <- data.frame(condition= sapply(colnames(DF), function(x) strsplit(x, "_")[[1]][1]), 
                            replicate= sapply(colnames(DF), function(x) strsplit(x, "_")[[1]][2]),
                            row.names = colnames(DF))
  dds <- DESeqDataSetFromMatrix(countData = DF,
                                colData = sampleTable,
                                design= as.formula("~ replicate + condition"))
  dds_object <- DESeq(dds)
  saveRDS(dds_object, "data/dds_object_1.0.rds")
}else
{
  dds_object <- readRDS("data/dds_object_1.0.rds")
}

# Compute fold Changes
comparisons <- CJ(c("W18", "W29", "WKD", "PH18", "PH29", "PHJ9", "PHJ11"), 
                  c("W18", "W29", "WKD", "PH18", "PH29", "PHJ9", "PHJ11"))
colnames(comparisons) <- c("nominator", "denominator")
comparisons <- comparisons[nominator!=denominator]
comparisons[, FC_file:= paste0("FC_tables/", nominator, "_vs_", denominator, "_FC.txt"), c(colnames(comparisons))]

diff <- comparisons[,
                    {
                      if(!file.exists(FC_file))
                      {
                        current <- as.data.table(as.data.frame(lfcShrink(dds_object, contrast= c("condition", nominator, denominator))), keep.rownames= T)
                        fwrite(current, FC_file, col.names=T, row.names = F, sep= "\t", quote= F, na= NA)
                      }
                      fread(FC_file)
                    }, c(colnames(comparisons))]

#--------------------------------------------------------------------#
# 3- MA plots
#--------------------------------------------------------------------#
diff[, exp:=  paste0(nominator, "_vs_", denominator), .(nominator, denominator)]
sub <- diff[exp %in% c("PH18_vs_W18", "PH29_vs_W29", "PHJ9_vs_WKD", "PHJ11_vs_WKD")]

pdf("pdf/MA_plots_PH_vs_W.pdf", width = 7, height = 8)
par(mfrow= c(2, 2), las= 1)
sub[, 
     {
       y <- log2FoldChange
       pch <- ifelse(abs(y)<=5, 16, 2)
       y[y >  5] <- 5
       y[y < -5] <- -5
       Cc <- ifelse(is.na(padj) | padj >= 0.01, "lightgrey", ifelse(y>0, "red", "blue"))
       plot(baseMean, y, col= Cc, ylim= c(-5, 5), pch= pch, cex= 0.5, log= "x", ylab= "log2FoldChange", xlab= "baseMean")
       mtext(exp[1], line = 1, cex= 0.7)
       leg <- c(paste(length(which(Cc=="red")), "Up"), paste(length(which(Cc=="blue")), "Down"), "(padj<0.01)")
       legend("topright", leg, bty= "n", text.col = c("red", "blue", "black"), cex = 0.8)
       print("DONE")
     }, exp]
dev.off()

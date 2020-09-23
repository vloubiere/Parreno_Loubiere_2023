setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
options(scipen= 99)
require(DESeq2)
require(gridExtra)
require(data.table)

# Import and clean counts
dat <- data.table(file= list.files("db/counts", ".txt", full.names = T, recursive = T))
dat[, sample:= gsub(".txt", "", basename(file))]
dat[, c("cdition", "rep"):= tstrsplit(sample, "_")]
dat[grepl("epiCancer", file), project:= "epiCancer"]
dat[grepl("SA2020/WRNAI|SA2020/PHRNAI", file), project:= "RNAI25"]
dat[is.na(project), project:= "AMM"]
dat <- dat[, fread(file), c(colnames(dat))]
colnames(dat)[6:7] <- c("symbol", "counts")
dat <- dat[!grepl("^__", symbol)]

#--------------------------------------------------------------------#
# 1- DESeq2 analysis
#--------------------------------------------------------------------#
# dds objects
dat[, 
    {
      current <- dcast(.SD, symbol~sample, value.var = "counts")
      dds_file <- paste0("Rdata/dds_object_", project, ".rds")
      if(!file.exists(dds_file))
      {
        DF <- data.frame(current[, -1], row.names = current$symbol)
        DF <- DF[rowSums(DF)>5,]
        sampleTable <- data.frame(condition= sapply(colnames(DF), function(x) strsplit(x, "_")[[1]][1]), 
                                  replicate= sapply(colnames(DF), function(x) strsplit(x, "_")[[1]][2]),
                                  row.names = colnames(DF))
        dds <- DESeqDataSetFromMatrix(countData = DF,
                                      colData = sampleTable,
                                      design= ~ replicate + condition)
        dds_object <- DESeq(dds)
        saveRDS(dds_object, dds_file)
      }
      dds_object <- readRDS(dds_file)
      
      # FC tables
      comparisons <- CJ(nominator= as.character(dds_object$condition), 
                        denominator= as.character(dds_object$condition), unique = T)
      comparisons <- comparisons[nominator!=denominator]
      comparisons <- comparisons[denominator %in% c("W18", "WKD", "W29", "PH18", "PH29", "WRNAI", "X2A", "X42D")]
      comparisons[, FC_file:= paste0("db/FC_tables/", nominator, "_vs_", denominator, "_FC.txt"), c(colnames(comparisons))]
      comparisons[, 
                  {
                    print(paste0(nominator, "x", denominator))
                    if(!file.exists(FC_file))
                    {
                      res <- lfcShrink(dds_object, contrast= c("condition", nominator, denominator))
                      current <- as.data.table(as.data.frame(res), keep.rownames= T)
                      fwrite(current, FC_file, col.names=T, row.names = F, sep= "\t", quote= F, na= NA)
                    }
                    print("DONE")
                  }, (comparisons)]
    }, project]

#--------------------------------------------------------------------#
# 2- MA plots
#--------------------------------------------------------------------#
diff <- data.table(file= list.files("db/FC_tables/", ".txt", full.names = T))
diff[, exp := gsub("_FC.txt", "", basename(file))]
diff <- diff[exp %in% c("PCXT1092A_vs_X2A", "PSCSUZ21B842D_vs_X42D", "EZ7312A_vs_X2A", "SUZ1212A_vs_X2A",
                        "PH18_vs_PH29", "PHJ9_vs_PH29", "PHJ11_vs_PH29",
                        "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18",
                        "PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29",
                        "PHRNAI_vs_WRNAI")]
diff <- diff[, fread(file), c(colnames(diff))]

pdf("pdf/MA_plots_transcriptomes.pdf", width = 7, height = 8)
par(mfrow= c(2, 2), las= 1)
diff[,
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

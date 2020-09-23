setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
options(scipen= 99)
require(DESeq2)
require(gridExtra)
require(data.table)

# Import and clean counts
dat <- fread("db/metadata/metadata_all.txt")
dat <- dat[method=="RNA"]
dat[!grepl("rep", rep), rep:= paste0("rep", rep)]
dat[, layout:= ifelse(.N==2, "PAIRED", "SINGLE"), id]
dat[project=="epiCancer" | grepl("RNAI$", sample), length:= 150]
dat[grepl("^WTE1416", sample), length:= 75]
dat[is.na(length), length:= 50]

# #--------------------------------------------------------------------#
# # 1- DESeq2 analysis
# #--------------------------------------------------------------------#
# current <- dat[, .(layout, length, rep, sample, id, counts_prefix)]
# current[, file:= list.files("db/counts/", paste0("^", basename(counts_prefix)), recursive = T, full.names = T), counts_prefix]
# 
# current <- current[, fread(file, col.names = c("FBgn", "counts")), current][!grepl("^__", FBgn)]
# dmat <- dcast(current, FBgn~id, value.var = "counts")
# DF <- data.frame(dmat[, -1], row.names = dmat$FBgn)
# DF <- DF[rowSums(DF)>5,]
# 
# sampleTable <- unique(dat[, .(id, sample, layout, length= factor(length), rep)])
# sampleTable[!id %in% colnames(DF), id:= grep(paste0("^X", id, "$"), colnames(DF), value = T), id]
# sampleTable <- data.frame(sampleTable[, -1], row.names = sampleTable$id)
# sampleTable <- sampleTable[match(colnames(DF), rownames(sampleTable)),]
# 
# dds <- DESeqDataSetFromMatrix(countData = DF,
#                               colData = sampleTable,
#                               design= ~ rep + sample)
#                               # design= ~ layout + length + rep + sample)
# dds_object <- DESeq(dds)
# saveRDS(dds_object, "Rdata/dds_object_all_transcriptomes.rds")


# FC tables
dds_object <- readRDS("Rdata/dds_object_all_transcriptomes.rds")
comparisons <- CJ(nominator= as.character(dds_object$sample), 
                  denominator= as.character(dds_object$sample), unique = T)
comparisons <- comparisons[nominator!=denominator]
comparisons <- comparisons[sapply(comparisons$nominator, function(x) dat$project[grep(x, dat$id)[1]]) == sapply(comparisons$denominator, function(x) dat$project[grep(x, dat$id)[1]])]
comparisons <- comparisons[!nominator %in% c("2A", "42D") & 
                           denominator %in% c("120hED", "WTE1416", "W18", "WKD", "W29", "PH18", "PH29", "WRNAI", "X2A", "X42D")]
comparisons[, FC_file:= paste0("db/FC_tables_all_transcriptomes/", nominator, "_vs_", denominator, "_FC.txt"), (comparisons)]
comparisons[, 
            {
              print(paste0(nominator, "x", denominator))
              if(!file.exists(FC_file))
              {
                res <- lfcShrink(dds_object, contrast= c("sample", nominator, denominator))
                current <- as.data.table(as.data.frame(res), keep.rownames= T)
                fwrite(current, FC_file, col.names=T, row.names = F, sep= "\t", quote= F, na= NA)
              }
              print("DONE")
            }, (comparisons)]

# #--------------------------------------------------------------------#
# # 2- MA plots
# #--------------------------------------------------------------------#
# diff <- data.table(file= list.files("db/FC_tables/", ".txt", full.names = T))
# diff[, exp := gsub("_FC.txt", "", basename(file))]
# diff <- diff[exp %in% c("PCXT1092A_vs_X2A", "PSCSUZ21B842D_vs_X42D", "EZ7312A_vs_X2A", "SUZ1212A_vs_X2A",
#                         "PH18_vs_PH29", "PHJ9_vs_PH29", "PHJ11_vs_PH29",
#                         "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18",
#                         "PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29",
#                         "PHRNAI_vs_WRNAI")]
# diff <- diff[, fread(file), c(colnames(diff))]
# 
# pdf("pdf/MA_plots_transcriptomes.pdf", width = 7, height = 8)
# par(mfrow= c(2, 2), las= 1)
# diff[,
#      {
#        y <- log2FoldChange
#        pch <- ifelse(abs(y)<=5, 16, 2)
#        y[y >  5] <- 5
#        y[y < -5] <- -5
#        Cc <- ifelse(is.na(padj) | padj >= 0.01, "lightgrey", ifelse(y>0, "red", "blue"))
#        plot(baseMean, y, col= Cc, ylim= c(-5, 5), pch= pch, cex= 0.5, log= "x", ylab= "log2FoldChange", xlab= "baseMean")
#        mtext(exp[1], line = 1, cex= 0.7)
#        leg <- c(paste(length(which(Cc=="red")), "Up"), paste(length(which(Cc=="blue")), "Down"), "(padj<0.01)")
#        legend("topright", leg, bty= "n", text.col = c("red", "blue", "black"), cex = 0.8)
#        print("DONE")
#      }, exp]
# dev.off()

# Compute WT activity
wt <- data.table(file= list.files("db/counts/RNA_development/", 
                                  "120hED", 
                                  full.names = T))
FPKM <- wt[, {
  .c <- readRDS(wt$file[1])
  res <- as.data.table(.c$counts, keep.rownames = "FBgn")
  colnames(res)[2] <- "FPKM"
  res$log2FPKM <- log2((res$FPKM+1)/(.c$annotation$Length/1000)/(sum(res$FPKM)/10e6))
  res
}, file]

FPKM <- FPKM[, .(log2FPKM= mean(log2FPKM)), FBgn]

# linear model to assess the activity of the motif
sel <- readRDS("Rdata/selected_motifs.rds")
counts <- readRDS("Rdata/all_motif_counts.rds")
counts <- counts[, colnames(counts) %in% c("FBgn", "symbol", unique(sel$motif)), with= F]
counts <- melt(counts, id.vars = c("FBgn", "symbol"))
counts <- counts[, .(value= log2(sum(value)+1)), .(FBgn, symbol, variable)]
counts[FPKM, log2FPKM:= i.log2FPKM, on= "FBgn"]
counts <- counts[log2FPKM>4]

act <- counts[, {
  .mod <- lm(log2FPKM~value)
  .sum <- summary(.mod)
  est <- summary(.mod)$coefficients[2,1]
  rsq <- summary(.mod)$r.squared
  pval <- pf(.sum$fstatistic[1],
             .sum$fstatistic[2],
             .sum$fstatistic[3],
             lower.tail=FALSE)
  data.table(estimate= est,
             rsquared= rsq,
             Fpval= pval)
},variable]

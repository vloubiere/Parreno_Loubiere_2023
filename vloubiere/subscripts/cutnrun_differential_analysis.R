setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")


# Import peaks files
K27Ac <- fread("db/peaks/K27Ac_cutnrun.bed", col.names = c("seqnames", "start", "end"), select = 1:3)
K27me3 <- fread("db/peaks/K27me3_cutnrun.bed", col.names = c("seqnames", "start", "end"), select = 1:3)
regions <- rbind(data.table(K27Ac, region_type= "K27Ac"),
                 data.table(K27me3, region_type= "K27me3"))

# Quantif signal
dat <- data.table(file= list.files("db/bw/cutnrun_merge/", full.names = T))
dat[, cdition:= gsub("_merge.bw", "", basename(file))]
dat <- dat[, cbind(regions, vl_bw_coverage(regions, file)), .(cdition, file)]

# dcast
dat <- dcast(dat, 
             seqnames+start+end+region_type~cdition, 
             value.var = "score")
dat[region_type=="K27me3", mean_act:= apply(.SD, 1, mean), .SDcols= patterns("H3K27me3")]
dat[region_type=="K27Ac", mean_act:= apply(.SD, 1, mean), .SDcols= patterns("H3K27Ac")]

# Compute Log2FC
cols <- grep("H3K27me3", names(dat), value = T)
dat[, paste0(cols, "_log2FC"):= lapply(.SD, function(x) c(log2(scale(x, center = F))-log2(scale(H3K27me3_PH18, center = F)))), .SDcols= cols]
cols <- grep("H3K27Ac", names(dat), value = T)
dat[, paste0(cols, "_log2FC"):= lapply(.SD, function(x) c(log2(scale(x, center = F))-log2(scale(H3K27Ac_PH18, center = F)))), .SDcols= cols]
cols <- c("seqnames", "start", "end", "region_type", "mean_act", grep("_log2FC$", names(dat), value= T))
dat <- dat[, ..cols]
dat$H3K27me3_PH18_log2FC <- NULL
dat$H3K27Ac_PH18_log2FC <- NULL
fwrite(dat, 
       "Rdata/K27_cutnrun_FC.txt")

# Hox?
dat[seqnames=="chr3R" & start<=16977663 & end>=16657408, HOX:= T]
dat[seqnames=="chr3R" & start<=7060439 & end>=6660880, HOX:= T]
dat[is.na(HOX), HOX:= F]

# Format for plotting MA plots
pl <- melt(dat, id.vars = c("seqnames", "start", "end", "region_type", "mean_act", "HOX"))
pl <- pl[(region_type=="K27me3" & grepl("K27me3", variable))
         | (region_type=="K27Ac" & grepl("K27Ac", variable))]
setorderv(pl, "HOX")

# Plot
pdf("pdf/cutnrun/MA_plots_K27_changes.pdf", width = 9, height = 6.75)
par(mfrow= c(2,3))
pl[, {
  name <- strsplit(as.character(variable[1]), "_")[[1]][2]
  pch <- ifelse(HOX, 21, 16)
  col <- adjustcolor(ifelse(abs(value)>1, "tomato", "lightgrey"), 0.5)
  col[HOX] <- "black"
  bg <- adjustcolor(ifelse(abs(value)>1, "tomato", "lightgrey"), 0.5)
  plot(log2(mean_act), 
       value,
       ylim= if(region_type=="H3K27me3") c(-3, 3) else c(-4,4),
       pch= pch,
       col= col,
       bg= bg,
       xlab= paste0("log2(", region_type, ")"),
       ylab= paste("log2FoldChange",  name, "vsv PH18"),
       main= variable[1],
       las= 1)
  legend("topright", 
         bty="n", 
         legend = paste("n=", length(which(value>1))))
  legend("bottomright", 
         bty="n", 
         legend = paste("n=", length(which(value<(-1)))))
  legend("right",
         0,
         bty="n",
         legend = paste("n=", length(which(between(value, -1, 1, incbounds = T)))))
  legend("topleft", 
         bty= "n", 
         pch= 21,
         legend = "HOX clusters")
  abline(h= -1, lty= 2)
  abline(h= 1, lty= 2)
}, .(region_type, variable)]
dev.off()


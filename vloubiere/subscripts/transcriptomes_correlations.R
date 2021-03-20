# Import
dat <- data.table(file= list.files("db/counts/", "W18|WKD|W29|Ez18|EzJ9|EzJ11|Ez29|PH18|PHJ9|PHJ11|PH29", full.names = T))
dat <- dat[, {.c <- as.data.table(readRDS(file)$counts, keep.rownames = "FBgn"); colnames(.c)[2] <- c("counts"); .c},  (dat)]
dat[, cdition:= gsub("_counts.rds", "", basename(file))]
dat <- dcast(dat, FBgn~cdition, value.var = "counts")

pdf("pdf/PCC.pdf", width = 10, height = 10)
pcc <- cor(as.matrix(dat, 1))
my_heatmap(pcc)
dev.off()

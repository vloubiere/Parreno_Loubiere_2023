# ALL
all <- data.table(file= list.files("db/counts/", recursive = T, full.names = T))
all <- all[, {
  .c <- as.data.table(readRDS(file)$counts, keep.rownames = "FBgn")
  colnames(.c)[2] <- c("counts")
  .c
}, file]

all[, cdition:= gsub("_counts.rds", "", basename(file))]
all <- dcast(all, FBgn~cdition, value.var = "counts")
pcc <- cor(as.matrix(all, 1))

pdf("pdf/PCC_all_epi_cancer.pdf", width = 20, height = 20)
par(mar= c(10,10,5,10), xaxs= "i", yaxs= "i")
vl_heatmap(pcc, display_numbers = T, col = matlab.like(10))
dev.off()

# EpiCancer
dat <- data.table(file= list.files("db/counts/RNA_epiCancer/", full.names = T))
dat <- dat[, {
  .c <- as.data.table(readRDS(file)$counts, keep.rownames = "FBgn")
  colnames(.c)[2] <- c("counts")
  .c}, (dat)]

dat[, cdition:= gsub("_counts.rds", "", basename(file))]
dat <- dcast(dat, FBgn~cdition, value.var = "counts")
pcc <- cor(as.matrix(dat, 1))

pdf("pdf/PCC_all_epi_cancer.pdf", width = 12, height = 12)
par(mar= c(10,10,5,10), xaxs= "i", yaxs= "i")
vl_heatmap(pcc, display_numbers = T)
dev.off()

pcc <- pcc[,!grepl("_T", colnames(pcc))]
pcc <- pcc[!grepl("_T", rownames(pcc)),]

pdf("pdf/PCC_pre_transplant_epi_cancer.pdf", width = 12, height = 12)
par(mar= c(10,10,5,10), xaxs= "i", yaxs= "i")
vl_heatmap(pcc, display_numbers = T)
dev.off()



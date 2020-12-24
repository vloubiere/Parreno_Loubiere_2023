setwd("/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
options(scipen = 99)
require(Rsubread)
require(DESeq2)
require(data.table)
require(kohonen)
require(plyr)
require(circlize)
require(rtracklayer)

if(!file.exits("Rdata/som_train_datasset.rds")){
  #-----------------------------------#
  # HiC
  #-----------------------------------#
  if(!file.exists("../../public_data/dm6/HiC/Rdata/HiC_score_all_chr.rds")){
    dat5 <- list.files("../../public_data/dm6/HiC/SA2020/", "5kb_10Mb", full.names = T)
    dat5 <- rbindlist(lapply(dat5, function(x) get(load(x))))
    dat5[, intervalID:= paste0(chrom1, "_", start1, "_", start2)]
    dat5 <- na.omit(dat5)
    dat10 <- list.files("../../public_data/dm6/HiC/SA2020/", "10kb_max", full.names = T)
    dat10 <- rbindlist(lapply(dat10, function(x) get(load(x))))
    dat10 <- rbind(dat10[, .(chrom1, start1, end1= round_any(end1-1, 5000 , floor), chrom2, start2, end2= round_any(end2-1, 5000 , floor), ED_HiC_score, LE_HiC_score)],
                   dat10[, .(chrom1, start1, end1= round_any(end1-1, 5000 , floor), chrom2, start2= round_any(start2+1, 5000, ceiling), end2, ED_HiC_score, LE_HiC_score)],
                   dat10[, .(chrom1, start1= round_any(start1+1, 5000, ceiling), chrom2, end1, start2, end2= round_any(end2-1, 5000 , floor), ED_HiC_score, LE_HiC_score)],
                   dat10[, .(chrom1, start1= round_any(start1+1, 5000, ceiling), chrom2, end1, start2= round_any(start2+1, 5000, ceiling), end2, ED_HiC_score, LE_HiC_score)])
    dat10[, intervalID:= paste0(chrom1, "_", start1, "_", start2)]
    dat <- rbindlist(list(dat5, dat10[!intervalID %in% dat5$intervalID]), idcol = T)
    dat[, norm:= 1-(ED_HiC_score/200+0.5)]
    dat[, coor1:= paste0(chrom1, ":", formatC(start1, width = 8, flag = "0", format = "d"), "-", formatC(end1, width = 8, flag = "0", format = "d"))]
    dat[, coor2:= paste0(chrom2, ":", formatC(start2, width = 8, flag = "0", format = "d"), "-", formatC(end2, width = 8, flag = "0", format = "d"))]
    dat <- dat[end1>start1 & end2>start2]
    saveRDS(dat, "../../public_data/dm6/HiC/Rdata/HiC_score_all_chr.rds")
  }
  if(!file.exists("../../public_data/dm6/HiC/Rdata/HiC_matrix.rds")){
    dat <- readRDS("../../public_data/dm6/HiC/Rdata/HiC_score_all_chr.rds")
    collapse <- rbindlist(mclapply(split(dat, dat$chrom1), function(x) {
      .c <- dcast(as.data.table(x), coor1~coor2, value.var = "norm")
      colnames(.c) <- c("coor", seq(ncol(.c)-1))
      .c <- as.matrix(.c, 1)
      .c[lower.tri(.c)] <- t(.c)[lower.tri(.c)]
      .c <- as.data.table(.c, keep.rownames= "coor")
      return(.c)
      }), fill= T)
    HiC <- as.matrix(collapse, 1)
    saveRDS(HiC, "../../public_data/dm6/HiC/Rdata/HiC_matrix.rds")
    # dat <- dat[start1>=6e6 & end1<=17.5e6 & start2>=6e6 & end2<=17.5e6]
    # bxd <- GRanges("chr3R", IRanges(16.25e6,17e6))
  }
  if(!exists("HiC")){
    HiC <- readRDS("../../public_data/dm6/HiC/Rdata/HiC_matrix.rds")
  }
  
  #-----------------------------------#
  # ChIP-Seq
  #-----------------------------------#
  ChIP <- as.data.table(GRanges(rownames(HiC)))
  bw <- list.files("../../public_data/dm6/bw/", "ED", full.names = T)
  bw <- bw[!grepl("FAIRE|INPUT", bw)]
  cols <- gsub("_merge|_ED|.bw", "", basename(bw))
  ChIP[, (cols):= mclapply(bw, function(x) {
    .c <- as.data.table(import.bw(x))
    .c <- .c[ChIP, sum(score*width)/sum(width), .EACHI, on= c("start<end", "end>start")]$V1
    .c[is.na(.c)] <- 0
    log2(.c+1)
  })]
  
  #-----------------------------------#
  # RNA-Seq
  #-----------------------------------#
  annot <- as.data.table(GRanges(rownames(HiC)))
  annot <- as.data.frame(annot[, .(GeneID= rownames(HiC), Chr= seqnames, Start= start, End= end, Strand= "+")])
  bam <- list.files("db/bam/", "^PH.*.bam$|^W.*.bam$", full.names = T)
  if(!file.exists("db/som/feat_coutns_HiC_bins.rds")){
    fc <- featureCounts(files = bam, annot.ext= annot, isGTFAnnotationFile = F, isPairedEnd = T, nthreads= 8)
    saveRDS(fc, "db/som/feat_coutns_HiC_bins.rds")
  }
  # Log2FoldChange
  if(!exists("FC")){
    FC <- readRDS("db/som/feat_coutns_HiC_bins.rds")
    FC <- as.data.table(FC$counts, keep.rownames = T)
    FC <- melt(FC, id.vars = "rn")
    FC[, c("cdition", "rep"):= tstrsplit(variable, "_|[.]", keep = c(1,2))]
    FC <- cbind(FC[, .(log2FC_PH18_vs_W18= mean(.SD[, log2(value[cdition=="PH18"]+1)-log2(value[cdition=="W18"]+1), rep]$V1)), rn],
                FC[, .(log2FC_PHJ9_vs_WKD= mean(.SD[, log2(value[cdition=="PHJ9"]+1)-log2(value[cdition=="WKD"]+1), rep]$V1)), rn],
                FC[, .(log2FC_PHJ11_vs_WKD= mean(.SD[, log2(value[cdition=="PHJ11"]+1)-log2(value[cdition=="WKD"]+1), rep]$V1)), rn],
                FC[, .(log2FC_PH29_vs_W29= mean(.SD[, log2(value[cdition=="PH29"]+1)-log2(value[cdition=="W29"]+1), rep]$V1)), rn])
    FC <- FC[, .SD, .SDcols = unique(names(FC))]
  }
  train <- list(chrom= factor(unlist(tstrsplit(rownames(train$HiC), ":", keep=1))),
                HiC= HiC, 
                ChIP= as.matrix(ChIP[, ATAC_Davie:SUZ12], rownames(HiC)),
                FC= as.matrix(FC, 1))
  saveRDS(train, "Rdata/som_train_datasset.rds")
}


#-----------------------------------#
# SOM
#-----------------------------------#
if(!file.exists("Rdata/som_test.rds")){
  train <- readRDS("Rdata/som_train_datasset.rds")
  grid <- somgrid(45, 45, "hexagonal", toroidal= F)
  som <- supersom(train, whatmap = 1:3, grid= grid, maxNA.fraction = 0.95)
  saveRDS(som, "Rdata/som_test.rds")
  }
if(!exists("som")){
  som <- readRDS("Rdata/som_test.rds")
}

#-----------------------------------#
# PLOT
#-----------------------------------#
Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))
coor <- as.data.table(GRanges(rownames(as.data.frame(som$data))))
coor[, chrom:= .GRP, seqnames]
pl <- cbind(som$data$ChIP, som$data$FC,
            coor[, .(chrom, start)], keep.rownames = T)
pl <- pl[, .SD, .SDcols = unique(names(pl))]
pl[, cl:= som$unit.classif]
pl <- melt(pl, id.vars = c("rn", "cl"))
pl <- pl[, .(value= mean(value)), keyby= .(variable, cl)]
empty <- seq(som$grid$xdim*som$grid$ydim)
empty <- empty[!empty %in% pl$cl]
empty <- data.table(variable= rep(unique(pl$variable), each= length(empty)), 
                    cl= rep(empty, length(unique(pl$variable))), value= NA)
pl <- rbind(pl, empty)
setorderv(pl, c("variable", "cl"))
  
pdf("pdf/som_plots.pdf", 17, 15)
par(mfrow=c(4,4))
plot(som, "counts", shape= "straight", border= NA)
pl[, {
  lim <- c(quantile(value, 0.025, na.rm= T), quantile(value, 0.975, na.rm= T))
  value[value<lim[1]] <- lim[1]
  value[value>lim[2]] <- lim[2]
  plot(som, "property", property= value, palette= Cc, shape= "straight", border= NA, 
       main= variable, zlim= lim, heatkeywidth = .5)
}, variable]
dev.off()
file.show("pdf/som_plots.pdf")



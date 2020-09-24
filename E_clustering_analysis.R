setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(scipen= 5)
require(data.table)
require(org.Dm.eg.db)
require(kohonen)
require(circlize)

# Kohonen
if(!file.exists("Rdata/som_genes_all.rds"))
{
  # dat
  diff <- data.table(file= list.files("db/FC_tables_all_transcriptomes/", ".txt", full.names = T))
  diff[, exp := gsub("_FC.txt", "", basename(file))]
  diff <- diff[, fread(file), (diff)]
  diff <- dcast(diff, rn~exp, value.var = "log2FoldChange")
  ChIP <- readRDS("Rdata/ChIP_normalized_enrichment_tss_body_cenh.rds")
  dat <- merge(diff, ChIP, by.x= "rn", by.y="FBgn")
  dat <- cbind(dat[,1], apply(dat[, -1], 2, function(x) 
  {
    lim <- quantile(x, c(0.01, 0.99), na.rm= T)
    x[x<lim[1]] <- lim[1]
    x[x>lim[2]] <- lim[2]
    return(x)
  }))
  
  # RNA
  RNAi <- dat[, .(PHJ9_vs_PH18, PHJ11_vs_PH18, PH29_vs_PH18,
                  PH18_vs_W18, PHJ9_vs_WKD, PHJ11_vs_WKD, PH29_vs_W29, PHRNAI_vs_WRNAI)]
  mut <- dat[, .(PCXT1092A_vs_2A, PSCSUZ21B842D_vs_42D, SUZ1212A_vs_2A)]
  dev <- dat[, .(L1ED_vs_E1416= `72hED_vs_120hED`, L2ED_vs_E1416= `96hED_vs_WTE1416`, L3ED_vs_E1416=  `120hED_vs_WTE1416`)]
  
  # ChIP
  prom <- dat[, .(tss_enr_POLIIGSE112868, tss_enr_PC, tss_enr_PH, tss_enr_PSC, tss_enr_SUZ12, tss_enr_EYGSE112868, 
                  tss_enr_H3K4me1, tss_enr_H3K4me2, tss_enr_H3K4me3, tss_enr_H4K20me1,
                  tss_enr_H3K27ac, tss_enr_H3K27me1, tss_enr_H3K27me2, tss_enr_H3K27me3, tss_enr_H2AK118Ub)]
  body <- dat[, .(gene_enr_H3K4me1, gene_enr_H4K20me1, gene_enr_H3K27ac, gene_enr_H3K27me1, gene_enr_H3K27me2, 
                  gene_enr_H3K27me3, gene_enr_H2AK118Ub)]
  enh <- dat[, .(enh_enr_POLIIGSE112868, enh_enr_PC, enh_enr_PH, enh_enr_PSC, enh_enr_SUZ12, enh_enr_EYGSE112868, 
                 enh_enr_H3K4me1, enh_enr_H3K4me2, enh_enr_H4K20me1,
                 enh_enr_H3K27ac, enh_enr_H3K27me1, enh_enr_H3K27me2, enh_enr_H3K27me3, enh_enr_H2AK118Ub)]
  
  # dat
  train <- list(RNAi, mut, dev, prom, body, enh)
  train <- lapply(train, function(x)
  {
    x <- as.matrix(x)
    rownames(x) <- dat$rn
    return(x)
  })
  
  mygrid <- somgrid(xdim= 28, ydim= 28, topo = 'hexagonal', toroidal = T)
  set.seed(1234)
  som.model <- supersom(train, user.weights = c(3,1,1,1,1,1), grid = mygrid, rlen = 500, maxNA.fraction = 0.3)
  saveRDS(som.model, "Rdata/som_genes_all.rds")
}

som <- readRDS("Rdata/som_genes_all.rds")
codes <- as.data.table(som$codes)
sel <- codes[, PHJ9_vs_PH18:PH29_vs_W29]
sub <- apply(sel, 1, function(x) any(abs(x)>1))
set.seed(1)
sel[(sub), kcl:= kmeans(.SD, 5)$cluster]
sel[is.na(kcl), kcl:= 6]

dat <- as.data.table(som$data)
dat[, cl:= som$unit.classif]
dat <- dat[, lapply(.SD, mean, na.rm= T), keyby= cl]
dat <- melt(dat, id.vars = "cl")
Cc <- colorRampPalette(c("cornflowerblue", "white", "white", "tomato"))


pdf("pdf/som_genes_all.pdf", 40, 10)
par(mfrow= c(4, 12))
dat[, 
    {
      lim <- c(-1, 1)
      value[value < lim[1]] <- lim[1]
      value[value > lim[2]] <- lim[2]
      plot(som, "property", property= value, palette= Cc, shape= "straight", zlim= lim, border= NA, main= variable)
      add.cluster.boundaries(som, sel$kcl, lwd= 1)
    }, variable]
dev.off()











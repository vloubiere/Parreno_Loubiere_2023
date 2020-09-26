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
  setkeyv(diff, "exp")
  diff <- diff[c("PCXT1092A_vs_2A", "PSCSUZ21B842D_vs_42D", "EZ7312A_vs_2A", "SUZ1212A_vs_2A",
                 "PH18_vs_PH29", "PHJ9_vs_PH29", "PHJ11_vs_PH29",
                 "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18",
                 "PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29", "PHRNAI_vs_WRNAI",
                 "WTE1416_vs_120hED", "72hED_vs_120hED", "96hED_vs_120hED",
                 "72hED_vs_WTE1416", "96hED_vs_72hED", "120hED_vs_96hED")]
  diff <- diff[, fread(file), (diff)]
  baseMean <- diff[exp=="PH29_vs_W29", .(rn, baseMean_PH29_vs_W29= log2(baseMean+1))]
  diff <- dcast(diff, rn~exp, value.var = "log2FoldChange")
  ChIP <- readRDS("Rdata/ChIP_normalized_enrichment_tss_body_cenh.rds")
  dat <- merge(diff, ChIP, by.x= "rn", by.y="FBgn")
  dat <- merge(dat, baseMean, by= "rn")
  dat <- melt(dat, id.vars = c("rn", "symbol"))
  dat[!is.na(value), value:= 
      {
        lim <- quantile(value, c(0.01, 0.99))
        value[value<lim[1]] <- lim[1]
        value[value>lim[2]] <- lim[2]
        value
      }, variable]
  dat <- dcast(dat, rn+symbol~variable, value.var = "value")
  
  # RNA
  RNAi <- dat[, .(PH18_vs_PH29, PHJ9_vs_PH29, PHJ11_vs_PH29,
                  PHJ9_vs_PH18, PHJ11_vs_PH18, PH29_vs_PH18,
                  PH18_vs_W18, PHJ9_vs_WKD, PHJ11_vs_WKD, PH29_vs_W29, PHRNAI_vs_WRNAI)]
  exp <- dat[, .(baseMean_PH29_vs_W29)]
  dev <- dat[, .(`WTE1416_vs_120hED`, `72hED_vs_120hED`, `96hED_vs_120hED`,
                 `72hED_vs_WTE1416`, `96hED_vs_72hED`, `120hED_vs_96hED`)]
  mut <- dat[, .(PCXT1092A_vs_2A, PSCSUZ21B842D_vs_42D, SUZ1212A_vs_2A)]
  
  # ChIP
  prom <- dat[, .(PC_tss, PH_tss, PSC_tss, SUZ12_tss, POLIIGSE112868_tss, EYGSE112868_tss, 
                  H3K4me1_tss, H3K4me2_tss, H3K4me3_tss, H4K20me1_tss, H3K36me3_tss,
                  H3K27ac_tss, H3K27me1_tss, H3K27me2_tss, H3K27me3_tss, H2AK118Ub_tss)]
  body <- dat[, .(H3K4me1_gene, H3K4me2_gene, H3K4me3_gene, H4K20me1_gene, H3K36me3_gene,
                  H3K27ac_gene, H3K27me1_gene, H3K27me2_gene, H3K27me3_gene, H2AK118Ub_gene)]
  enh <- dat[, .(PC_cenh, PH_cenh, PSC_cenh, SUZ12_cenh, POLIIGSE112868_cenh, EYGSE112868_cenh, 
                  H3K4me1_cenh, H3K4me2_cenh, H3K4me3_cenh, H4K20me1_cenh, H3K36me3_cenh,
                  H3K27ac_cenh, H3K27me1_cenh, H3K27me2_cenh, H3K27me3_cenh, H2AK118Ub_cenh)]
  # dat
  train <- list(RNAi, do.call(cbind, list(exp, dev, mut, prom, body, enh)))
  train <- lapply(train, function(x)
  {
    x <- as.matrix(x)
    rownames(x) <- dat[, paste0(rn, "__", symbol)]
    return(x)
  })
  
  mygrid <- somgrid(xdim= 28, ydim= 28, topo = 'hexagonal', toroidal = T)
  set.seed(1234)
  som.model <- supersom(train, user.weights = c(5,1), grid = mygrid, rlen = 500, maxNA.fraction = 0.6)
  saveRDS(som.model, "Rdata/som_genes_all.rds")
}

som <- readRDS("Rdata/som_genes_all.rds")
codes <- as.data.table(som$codes)
sel <- codes[, PH18_vs_PH29:PHRNAI_vs_WRNAI]
sub <- apply(sel, 1, function(x) any(abs(x)>1))
set.seed(1)
sel[(sub), kcl:= kmeans(.SD, 5)$cluster]
sel[is.na(kcl), kcl:= 6]

dat <- as.data.table(som$data)
dat[, cl:= som$unit.classif]
dat <- dat[, lapply(.SD, mean, na.rm= T), keyby= cl]
dat <- melt(dat, id.vars = "cl")

pdf("pdf/som_genes_all.pdf", 40, 10)
par(mfrow= c(4, 12))
dat[, 
    {
      lim <- quantile(value, c(0.01, 0.99), na.rm = T)
      if(lim[1]<0)
      {
        lim <- c(-max(abs(lim)), max(abs(lim)))
      }else
      {
        lim <- c(0, lim[2])
      }
      value[value<lim[1]] <- lim[1]
      value[value>lim[2]] <- lim[2]
      Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))
      plot(som, "property", property= value, palette= Cc, shape= "straight", zlim= lim, border= NA, main= variable)
      add.cluster.boundaries(som, sel$kcl, lwd= 1)
    }, variable]
dev.off()


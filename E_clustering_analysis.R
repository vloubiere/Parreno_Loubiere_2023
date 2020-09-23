setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(scipen= 5)
require(data.table)
require(org.Dm.eg.db)
require(kohonen)
require(circlize)

# Kohonen
if(!file.exists("Rdata/som_genes.rds"))
{
  diff <- data.table(file= list.files("db/FC_tables/", ".txt", full.names = T))
  diff[, exp := gsub("_FC.txt", "", basename(file))]
  diff <- diff[exp %in% c("PCXT1092A_vs_X2A", "PSCSUZ21B842D_vs_X42D", "EZ7312A_vs_X2A", "SUZ1212A_vs_X2A",
                          "PH18_vs_PH29", "PHJ9_vs_PH29", "PHJ11_vs_PH29",
                          "PHJ9_vs_PH18", "PHJ11_vs_PH18", "PH29_vs_PH18",
                          "PH18_vs_W18", "PHJ9_vs_WKD", "PHJ11_vs_WKD", "PH29_vs_W29",
                          "PHRNAI_vs_WRNAI")]
  diff <- diff[, fread(file), (diff)]
  diff <- dcast(diff, rn~exp, value.var = "log2FoldChange")
  
  ChIP <- readRDS("Rdata/ChIP_normalized_enrichment_tss_body.rds")
  
  dat <- merge(diff, ChIP, by.x= "rn", by.y= "FBgn")
  mat <- as.matrix(dat, 1)
  mat <- apply(mat, 2, function(x) 
  {
    lim <- quantile(x, c(0.01, 0.99), na.rm= T)
    x[x<lim[1]] <- lim[1]
    x[x>lim[2]] <- lim[2]
    return(x)
  })
  
  mygrid <- somgrid(xdim= 28, ydim= 28, topo = 'hexagonal', toroidal = T, )
  set.seed(1234)
  som.model <- supersom(mat, grid = mygrid, rlen = 500)
  saveRDS(som.model, "Rdata/som_genes.rds")
}

som <- readRDS("Rdata/som_genes.rds")
# codes <- as.data.table(som$codes)
# dat <- as.data.table(som$data, keep.rownames = T)
dat <- as.data.table(som$data)
dat[, cl:= som$unit.classif]
dat <- dat[, lapply(.SD, mean), keyby= cl]
dat <- melt(dat, id.vars = "cl")
Cc <- colorRampPalette(c("cornflowerblue", "white", "white", "tomato"))

pdf("pdf/som_genes.pdf", 40, 10)
par(mfrow= c(4, 12))
dat[, 
    {
      lim <- max(abs(range(value)))
      plot(som, "property", property= value, palette= Cc, shape= "straight", zlim= c(-lim, lim), border= NA, main= variable)
    }, variable]
dev.off()



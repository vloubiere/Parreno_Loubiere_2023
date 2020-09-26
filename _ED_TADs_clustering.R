setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
require(pheatmap)
require(rtracklayer)
require(data.table)
require(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

# Import TADs
tad <- fread("Rdata/ED_TADS_2kb.txt")
colnames(tad)[4] <- "name"
tad <- GRanges(tad[order(chrom, start)])

# ChIP data 
ChIP <- data.table(file= list.files("db/bed/ChIP", ".bed", full.names = T))
ChIP[, cdition:= gsub("_uniq.bed|_rep1|_rep2", "", basename(file))]
ChIP[cdition %in% c("PC", "PH", "H3K27me3"), input_group:= "INPUTa"]
ChIP[cdition %in% c("PSC", "SUZ12"), input_group:= "INPUTv"]
ChIP[cdition %in% c("EYGSE112868", "POLIIGSE112868"), input_group:= "INPUTGSE112868"]
ChIP[grepl("^H3|^H4|^H2", cdition) & is.na(input_group), input_group:= "INPUTh"]
if(!file.exists("Rdata/TAD_ChIP_quantif.rds"))
{
  quantif <- ChIP[, my_countReads(tad, file), (ChIP)]
  saveRDS(quantif, "Rdata/TAD_ChIP_quantif.rds")
}
if(!file.exists("Rdata/som_TADs_ED.rds"))
{
  quantif <- readRDS("Rdata/TAD_ChIP_quantif.rds")
  norm <- quantif[, .(norm_counts= (sum(counts)+1)/sum(total_reads)*1e6), .(name, cdition, input_group)]
  norm[norm, input_norm_counts:= i.norm_counts, on= c("name", "input_group==cdition")]
  norm <- norm[!grepl("INPUT", cdition)]
  norm[, enr:= log2(norm_counts)-log2(input_norm_counts)]
  dmat <- dcast(norm, name~cdition, value.var = "enr")
  setcolorder(dmat, c("name", "PC", "PH", "PSC", "SUZ12", "POLIIGSE112868", "EYGSE112868", "H3K4me3", "H3K4me2", 
                      "H3K36me3", "H4K20me1", "H3K27ac", "H3K4me1", "H3K27me1", "H3K27me3", "H2AK118Ub", "H3K27me2"))
  colnames(dmat)[6] <- "PolII"
  colnames(dmat)[7] <- "EY"
  
  mygrid <- somgrid(xdim= 3, ydim= 3, topo = 'hexagonal', toroidal = T)
  mat <- scale(as.matrix(dmat, 1))
  set.seed(1234)
  som.model <- supersom(mat, grid = mygrid, rlen = 500, maxNA.fraction = 0.6)
  
  Cc <- colorRampPalette(c("black", "grey20", "grey40", "grey60", "gold", "darkgoldenrod3", "darkred", "red", "tomato", 
                           "seagreen4", "limegreen", "chartreuse3", "palegreen", "cornflowerblue", "royalblue2", "lightgrey"))
  plot(som.model, "codes", shape= "straight", codeRendering= "segments", palette= Cc)
  plot(som.model, "counts", shape= "straight")
  
  saveRDS(som.model, "Rdata/som_TADs_ED.rds")
}
setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(VennDiagram)

#-------------------#
dat <- readRDS("Rdata/final_FC_table.rds")
dat[padj<0.01 & log2FoldChange>0, class:= "up"]
dat[padj<0.01 & log2FoldChange<0, class:= "down"]
dat <- dat[!is.na(class)]
setkeyv(dat, c("cdition", "class"))

# PLOT
pdf("pdf/genes_overlaps/overlap_up_down_genes.pdf")
par(xaxs= "i",
    yaxs= "i",
    mar= rep(0, 4))
for(cdition in c("29_TS", "D11_TS", "D9_TS"))
{
  c1 <- paste0("PH", cdition)
  c2 <- paste0("EZ", cdition)
  
  for(FC in c("up", "down"))
  {
    plot.new()
    text(0.5, 0.95, paste0(c1, " VS ", c2, "  ", FC))
    VennDiagram::draw.pairwise.venn(nrow(dat[.(c1, FC)]), 
                                    nrow(dat[.(c2, FC)]), 
                                    length(which(dat[.(c1, FC), FBgn] %in% dat[.(c2, FC), FBgn])), 
                                    cex = 3,
                                    fill = c("#66CC99", "#FFFFCC"))
    
  }
}
dev.off()

#----------------#
cmb <- data.table(c1= c("PH29_vs_W29", "PH29_vs_W29", "PHD9_vs_WKD", "EZ29_vs_W29", "EZ29_vs_W29", "EZD9_vs_WKD"),
                  c2= c("PHD9_vs_WKD", "PHD11_vs_WKD", "PHD11_vs_WKD", "EZD9_vs_WKD", "EZD11_vs_WKD", "EZD11_vs_WKD"))
cmb <- cmb[, .(c3= c("up", "down")), cmb]

pdf("pdf/genes_overlaps/overlap_timepoints_genes.pdf")
par(xaxs= "i",
    yaxs= "i",
    mar= rep(0, 4))
cmb[, {
  plot.new()
  text(0.5, 0.95, paste0(c1, " VS ", c2, "  ", c3))
  VennDiagram::draw.pairwise.venn(nrow(dat[cdition==c1 & class==c3]), 
                                  nrow(dat[cdition==c2 & class==c3]), 
                                  length(which(dat[cdition==c1 & class==c3, V1] %in% dat[cdition==c2 & class==c3, V1])), 
                                  cex = 3)
  print("")
}, (cmb)]
dev.off()

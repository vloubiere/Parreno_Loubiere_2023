setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

#-------------------------#
# Import ph29 vs ph18 dose comparison
#-------------------------#
dat <- data.table(file= c("db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH18_vs_RNA_W18.txt",
                          "db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH29_vs_RNA_W29.txt"))
dat[, cdition:= c("PH18", "PH29")]
dat <- dat[, fread(file), (dat)]
#-------------------------#
# Add symbols
#-------------------------#
symbols <- unique(readRDS("Rdata/final_FC_table.rds")[, 2:3])
dat[symbols, symbol:= i.symbol, on= "V1==FBgn"]
dat <- dat[, .(FBgn= V1, symbol, log2FoldChange, padj, cdition)]
dat <- dcast(dat, 
             FBgn+symbol~cdition, 
             value.var = c("log2FoldChange", "padj"))
dat[, Cc1:= "lightgrey"]
dat[, Cc2:= "lightgrey"]
dat[log2FoldChange_PH18>1 & padj_PH18<=0.05 & padj_PH29>0.05, c("Cc1", "cl"):= .("tomato", 1)]
dat[log2FoldChange_PH18<(-1) & padj_PH18<=0.05 & padj_PH29>0.05, c("Cc1", "cl"):= .("cornflowerblue", 2)]
dat[log2FoldChange_PH29>1 & padj_PH29<=0.05 & padj_PH18>0.05, c("Cc2", "cl"):= .("gold", 3)]
dat[log2FoldChange_PH29<(-1) & padj_PH29<=0.05 & padj_PH18>0.05, c("Cc2", "cl"):= .("green", 4)]
dat[, Cc:= colorRampPalette(c(Cc1, Cc2))(3)[2], .(Cc1, Cc2)]

#-------------------------#
# Scatterplot
#-------------------------#
dir.create("pdf/comparison_ph18_ph29_DOSE", 
           showWarnings = F)
pdf("pdf/comparison_ph18_ph29_DOSE/scatterplot_plot_ph18_vs_29.pdf", 
    width = 4.5)
par(pty= "s")
plot(dat$log2FoldChange_PH18, 
     dat$log2FoldChange_PH29, 
     col= adjustcolor(dat$Cc, 0.6),
     pch= 19,
     cex= 0.6)
abline(0, 1)
dev.off()

#-------------------#
# GO
#-------------------#
pdf("pdf/comparison_ph18_ph29_DOSE/GO_ph18_vs_29.pdf", 
    width =  10)
.l <- split(dat$FBgn, dat$cl)
vl_GO_clusters(.l, 
               N_top = 10, 
               padj_cutoff = 0.05, 
               cex = 0.6)
dev.off()



#-------------------#
# Protein networks
#-------------------#

net <- vl_STRING_interaction(sub[, symbol], 
                            size = abs(sub[, log2FoldChange]), 
                            cex.label = abs(sub[, log2FoldChange]), 
                            col = ifelse(sub$log2FoldChange>0, "tomato", "cornflowerblue"),
                            score_cutoff = 700)

pdf("pdf/comparison_ph18_ph29_DOSE/protein_network_ph18_vs_29.pdf")
vl_STRING_network(net, 
                  cex.vertices = 2,
                  cex.vertices.labels = 0.1)
dev.off()

files <- list("db/FC_tables/RNA_epiCancer_Ez18_vs_W18_FC.txt",
            "db/FC_tables/RNA_epiCancer_EzJ9_vs_WKD_FC.txt",
            "db/FC_tables/RNA_epiCancer_EzJ11_vs_WKD_FC.txt",
            "db/FC_tables/RNA_epiCancer_Ez29_vs_W29_FC.txt",
            "db/FC_tables/RNA_epiCancer_PH18_vs_W18_FC.txt",
            "db/FC_tables/RNA_epiCancer_PHJ9_vs_WKD_FC.txt",
            "db/FC_tables/RNA_epiCancer_PHJ11_vs_WKD_FC.txt",
            "db/FC_tables/RNA_epiCancer_PH29_vs_W29_FC.txt",
            "db/FC_tables/RNA_Paro_2018_RNA_4WED_vs_RNA_WTED_FC.txt",
            "db/FC_tables/RNA_Paro_2018_RNA_8WED_vs_RNA_WTED_FC.txt",
            "db/FC_tables/RNA_Paro_2018_RNA_14WED_vs_RNA_WTED_FC.txt",
            "db/FC_tables/RNA_Paro_2018_RNA_ecdED_vs_RNA_WTED_FC.txt", 
            "db/FC_tables/RNA_phRNAi_SA2020_RNA_PHRNAI_ED_vs_RNA_WRNAI_ED_FC.txt", 
            "db/FC_tables/RNA_mutants_SA2020_RNA_EZ7312A_ED_vs_RNA_2A_ED_FC.txt", 
            "db/FC_tables/RNA_mutants_SA2020_RNA_PCXT1092A_ED_vs_RNA_2A_ED_FC.txt", 
            "db/FC_tables/RNA_mutants_SA2020_RNA_PSCSUZ21B842D_ED_vs_RNA_42D_ED_FC.txt",
            "db/FC_tables/RNA_mutants_SA2020_RNA_SUZ1212A_ED_vs_RNA_2A_ED_FC.txt",
            "external_data/dlg_transcriptome_Bilder_2015_elife-03189-supp2-v2.xlsx",
            "external_data/scrib_transcriptome_Bilder_2015_elife-03189-supp1-v2.xlsx",
            "external_data/PSC_transcriptome_Bilder_2015_elife-03189-supp3-v2.xlsx")
dat <- lapply(files, function(x) 
{
  if(grepl(".txt$", x))
    x <- fread(x)
  else if(grepl(".xlsx", x))
  {
    x <- as.data.table(read_xlsx(x))
    colnames(x)[7] <- "log2FoldChange"
  }
  colnames(x)[1] <- "FBgn"
  return(x)
})
names(dat) <- files
dat <- rbindlist(dat, 
                 idcol = "cdition", 
                 fill= T)
dat[, cdition:= gsub("_ED|ED|_FC.txt$", "", basename(cdition)), cdition]
dat[, cdition:= gsub("2018_|RNA_", "", cdition), cdition]
dat[, cdition:= gsub("phRNAi_SA2020", "25C", cdition), cdition]
dat[, cdition:= gsub("mutants_SA2020", "CpGmut", cdition), cdition]
dat[, cdition:= gsub("^(.*)_transcriptome_Bilder.*", "Bilder_\\1", cdition), cdition]
dat[, log2FoldChange:= as.numeric(log2FoldChange)]
dat[, c("min", "max"):= as.list(quantile(log2FoldChange[is.finite(log2FoldChange)], c(0.01, 0.99), na.rm = T)), cdition]
dat[log2FoldChange>max, log2FoldChange:= max]
dat[log2FoldChange<min, log2FoldChange:= min]
.m <- dcast(dat, FBgn~cdition, value.var = "log2FoldChange")
.m <- as.matrix(na.omit(.m), 1)
.m <- .m[, c("Bilder_PSC",
             "Bilder_dlg", 
             "Bilder_scrib", 
             "CpGmut_EZ7312A_vs_2A", 
             "CpGmut_PCXT1092A_vs_2A", 
             "CpGmut_PSCSUZ21B842D_vs_42D",
             "CpGmut_SUZ1212A_vs_2A", 
             "Paro_4W_vs_WT", 
             "Paro_8W_vs_WT", 
             "Paro_14W_vs_WT", 
             "Paro_ecd_vs_WT", 
             "epiCancer_Ez18_vs_W18", 
             "epiCancer_EzJ11_vs_WKD", 
             "epiCancer_EzJ9_vs_WKD", 
             "epiCancer_Ez29_vs_W29", 
             "epiCancer_PH18_vs_W18", 
             "epiCancer_PHJ11_vs_WKD", 
             "epiCancer_PHJ9_vs_WKD", 
             "epiCancer_PH29_vs_W29", 
             "25C_PHRNAI_vs_WRNAI")]
pca <- prcomp(scale(.m, scale = F))
res <- as.data.table(pca$rotation[, c("PC1", "PC2")], keep.rownames = T)
res[, Cc:= c("orchid1",
             "darkorchid1",
             "darkorchid4",
             "gold",
             "goldenrod",
             "goldenrod4",
             "sienna4",
             "olivedrab1",
             "limegreen",
             "olivedrab3",
             "olivedrab4",
             "lightsteelblue1",
             "cornflowerblue",
             "blue",
             "navy",
             "pink",
             "indianred2", 
             "red2", 
             "red3",
             "black")]

pdf("pdf/PCA_transcriptomes.pdf", width = 9, height = 5.5)
par(mar= c(4,4,2,21), 
    las= 1)
plot(res$PC1, 
     res$PC2, 
     col= res$Cc, 
     pch= 19, 
     xlab= "PC1",
     ylab= "PC2",
     cex= 1.5)
legend(par("usr")[2], 
       par("usr")[4]+0.025, 
       bty= "n", 
       col= res$Cc, 
       legend = res$rn, 
       pch= 19,
       xpd=T,
       cex= 1)
dev.off()

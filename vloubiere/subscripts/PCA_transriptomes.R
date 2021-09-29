setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(diagram)

dat <- data.table(file= c("db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH18_vs_RNA_W18.txt",
                          "db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH21_vs_RNA_W21.txt",
                          "db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH25_vs_RNA_W25.txt",
                          "db/FC_tables/RNA_epiCancer_dose_ED_handDissect_Parreno_RNA_PH29_vs_RNA_W29.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_EZ18_vs_RNA_W18.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_EZD11_vs_RNA_WKD.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_EZD9_vs_RNA_WKD.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_EZ29_vs_RNA_W29.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_PH18_vs_RNA_W18.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_PHD11_vs_RNA_WKD.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_PHD9_vs_RNA_WKD.txt",
                          "db/FC_tables/RNA_epiCancer_ED_handDissect_Parreno_RNA_PH29_vs_RNA_W29.txt",
                          "db/FC_tables/RNA_Paro_2018_NA_Paro_RNA_4WED_vs_RNA_WTED.txt",
                          "db/FC_tables/RNA_Paro_2018_NA_Paro_RNA_8WED_vs_RNA_WTED.txt",
                          "db/FC_tables/RNA_Paro_2018_NA_Paro_RNA_14WED_vs_RNA_WTED.txt",
                          "db/FC_tables/RNA_Paro_2018_NA_Paro_RNA_ecdED_vs_RNA_WTED.txt", 
                          "db/FC_tables/RNA_phRNAi_SA2020_ED_NA_Loubiere_RNA_PHRNAI_ED_vs_RNA_WRNAI_ED.txt", 
                          "db/FC_tables/RNA_mutants_SA2020_ED_NA_Martinez_RNA_EZ7312A_ED_vs_RNA_2A_ED.txt", 
                          "db/FC_tables/RNA_mutants_SA2020_ED_NA_Martinez_RNA_SUZ1212A_ED_vs_RNA_2A_ED.txt",
                          "db/FC_tables/RNA_mutants_SA2020_ED_NA_Martinez_RNA_PCXT1092A_ED_vs_RNA_2A_ED.txt", 
                          "db/FC_tables/RNA_mutants_SA2020_ED_NA_Martinez_RNA_PSCSUZ21B842D_ED_vs_RNA_42D_ED.txt"))
dat <- dat[, .(cdition= gsub(".*_RNA_(.*)_vs_.*", "\\1", basename(file))), file]
dat[grepl("RNA_epiCancer", file), cdition:= paste0(cdition, ifelse(grepl("_dose_", file), "_dose", "_TS")), file]

# Import, Format and Clip log2FC
dat <- dat[, fread(file), (dat)]
dat[, log2FoldChange:= as.numeric(log2FoldChange)]
dat <- dat[!is.na(log2FoldChange)]
dat[, c("min", "max"):= as.list(quantile(log2FoldChange[is.finite(log2FoldChange)], c(0.005, 0.995), na.rm = T)), cdition]
dat[log2FoldChange>max, log2FoldChange:= max]
dat[log2FoldChange<min, log2FoldChange:= min]

# Cast matrix
.m <- dcast(dat, V1~cdition, value.var = "log2FoldChange")
.m <- as.matrix(na.omit(.m), 1)

# PCA
pca <- prcomp(scale(.m))
res <- as.data.table(pca$rotation[, c("PC1", "PC2")], keep.rownames = T)
res <- res[match(rn, colnames(.m))]
res[, rn:= factor(rn, levels = unique(dat$cdition))]
setkeyv(res, "rn")
res["PH18_dose", Cc:= "orchid1"]
res["PH21_dose", Cc:= "darkorchid1"]
res["PH25_dose", Cc:= "purple"]
res["PH29_dose", Cc:= "darkorchid4"]
res["4WED", Cc:= "olivedrab1"]
res["8WED", Cc:= "limegreen"]
res["14WED", Cc:= "olivedrab3"]
res["ecdED", Cc:= "olivedrab4"]
res["EZ18_TS", Cc:= "lightsteelblue1"]
res["EZD11_TS", Cc:= "cornflowerblue"]
res["EZD9_TS", Cc:= "blue"]
res["EZ29_TS", Cc:= "navy"]
res["PH18_TS", Cc:= "pink"]
res["PHD11_TS", Cc:= "pink3"]
res["PHD9_TS", Cc:= "indianred2"]
res["PH29_TS", Cc:= "red2"]
res["PHRNAI_ED", Cc:= "black"]
res["EZ7312A_ED", Cc:= "gold"]
res["SUZ1212A_ED", Cc:= "goldenrod"]
res["PCXT1092A_ED", Cc:= "goldenrod4"]
res["PSCSUZ21B842D_ED", Cc:= "sienna4"]
res[grepl("TS$", rn), pch:= 17]
res[grepl("dose$", rn), pch:= 15]
res[is.na(pch), pch:= 19]

pdf("pdf/PCA_transcriptomes.pdf", width = 9.25, height = 5.5)
par(mar= c(4,4,2,21), 
    las= 1)
plot(res$PC1, 
     res$PC2, 
     col= res$Cc, 
     xlab= "PC1",
     ylab= "PC2",
     cex= 1.5, 
     xlim= c(-0.05, 0.37),
     pch= res$pch)
legend(par("usr")[2], 
       par("usr")[4]+0.025, 
       bty= "n", 
       col= res$Cc, 
       legend = res$rn, 
       pch= res$pch,
       xpd=T,
       cex= 1)
setkeyv(res, "rn")
dev.off()

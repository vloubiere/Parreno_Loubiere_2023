setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

# Import metadata
dat <- read_xlsx("Rdata/raw_metadata_final.xlsx")
dat <- as.data.table(dat)
dat <- dat[grepl("RNA_epiCancer", project) & tissue=="ED" & method=="handDissect"]
# Find count files
dat[, 
    counts_file:= paste0("db/counts/", project, "/", prefix_output, "_counts_only.txt"), .(project, prefix_output)]
# Build cdition names
dat[, group:= ifelse(grepl("_dose", counts_file), "dose", "TS")]
dat <- dat[, {
        .v <- c(gsub("RNA_", "", cdition), rep, group)
        .(pl= paste0(.v, collapse = "_"))
}, .(cdition, rep, group, counts_file)]
vCc <- matlab.like2(max(dat[group=="TS", .GRP, cdition]$GRP))
dat[group=="TS", Cc:= vCc[.GRP], cdition]
vCc <- matlab.like2(max(dat[group=="dose", .GRP, cdition]$GRP))
dat[group=="dose", Cc:= vCc[.GRP], cdition]
dat <- dat[, fread(counts_file), (dat)]
# Correlation matrix
mat <- dcast(dat, 
             FBgn~pl, 
             value.var = "counts")
mat <- as.matrix(mat, 1)
mat <- mat[rowSums(mat)>20,]
cor <- cor(mat)
# Points formating
leg <- unique(dat[colnames(mat), .(pl, Cc), on= "pl"])
leg[grepl("^W.*_TS$", pl), pch:= 21]
leg[grepl("^W.*_dose$", pl), pch:= 19]
leg[grepl("^PH.*_TS$", pl), pch:= 15]
leg[grepl("^PH.*_dose$", pl), pch:= 22]
leg[grepl("^EZ.*_TS$", pl), pch:= 17]
leg[grepl("^EZ.*_dose$", pl), pch:= 25]

# Heatmap
dir.create("pdf/PCC_replicates",
           showWarnings = F)
pdf("pdf/PCC_replicates/RNA_epiCancer_Parreno_heatmap.pdf",
    width = 9,
    height = 9)
par(cex.axis= 0.3)
vl_heatmap(cor)
dev.off()

# Dendrograms
dd <- as.dist(1 - cor(mat, method= "pearson"))
hc <- hclust(dd, method = "ward.D2")
pdf("pdf/PCC_replicates/RNA_epiCancer_Parreno_dendro.pdf", 
    width = 10, 
    height = 5)
par(bg= "grey")
plot(hc, 
     hang= -1, 
     cex= 0.6, 
     las= 1, 
     xlab= "Samples replicates",
     lend= 2)
points(seq(hc$order),
       rep(0, length(hc$order)),
       col= leg$Cc[hc$order],
       pch= leg$pch[hc$order],
       bg= "grey",
       cex= 0.6)
dev.off()

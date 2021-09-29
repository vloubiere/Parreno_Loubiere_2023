dat <- data.table(file= list.files("db/fpkms/", full.names = T))
dat <- dat[, {
  .c <- fread(file)
  melt(.c, id.vars = "FBgn")
}, file]
dat[, cdition:= gsub("_counts.rds", "", variable), variable]
dat <- dcast(dat, FBgn~cdition, value.var = "value")

mat <- as.matrix(dat, 1)
mat[is.na(mat)] <- 0
mat <- log2(mat+0.001)
# mat <- scale(mat)
pca <- prcomp(mat)
res <- as.data.table(pca$rotation[, c("PC1", "PC2")], keep.rownames= "cdition")
res[grepl("Ez18", cdition), Cc:= "lightsteelblue1"]
res[grepl("EzJ9", cdition), Cc:= "cornflowerblue"]
res[grepl("EzJ11", cdition), Cc:= "blue"]
res[grepl("Ez29", cdition), Cc:= "navy"]
res[grepl("PH18", cdition), Cc:= "pink3"] #pink
res[grepl("PHJ9", cdition), Cc:= "indianred2"]
res[grepl("PHJ11", cdition), Cc:= "red2"]
res[grepl("PH29", cdition), Cc:= "red3"] # darkred
res[grepl("E1416", cdition), Cc:= "orchid1"]
res[grepl("72hED", cdition), Cc:= "darkorchid1"]
res[grepl("96hED", cdition), Cc:= "darkorchid4"]
res[grepl("120hED", cdition), Cc:= "purple"]
res[grepl("WTED", cdition), Cc:= "olivedrab1"]
res[grepl("4WED", cdition), Cc:= "green"]
res[grepl("8WED", cdition), Cc:= "limegreen"]
res[grepl("14WED", cdition), Cc:= "olivedrab3"]
res[grepl("ecdED", cdition), Cc:= "olivedrab4"]
res[grepl("W18", cdition), Cc:= "grey50"]
res[grepl("W29", cdition), Cc:= "grey90"]
res[grepl("WKD", cdition), Cc:= "grey70"]
res[grepl("WRNAI", cdition), Cc:= "antiquewhite3"]
res[grepl("PHRNAI", cdition), Cc:= "black"]
res[grepl("_2A_", cdition), Cc:= "gold"]
res[grepl("PCXT", cdition), Cc:= "goldenrod"]
res[grepl("EZ731", cdition), Cc:= "goldenrod4"]
res[grepl("SUZ12", cdition), Cc:= "sienna4"]
res[grepl("PSCSUZ2", cdition), Cc:= "pink"]
res[grepl("_42D_", cdition), Cc:= "darkred"]


par(mar= c(5,5,2,10))
plot(res[, PC1], 
     res[, PC2], 
     col= res$Cc, 
     pch= 19)
leg <- unique(res[, .(Cc, gsub("(.*)_.*$", "\\1", cdition))])
legend(par("usr")[2], 
       par("usr")[4],
       bty= "n",
       pch= 19,
       col= leg$Cc, 
       legend = leg$V2,
       xpd= T)

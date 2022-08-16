require(data.table)
require(flowCore)
require(sp)

dat <- data.table(flowCore::read.FCS("db/fcs/D 2017-06-15_15062017 PSCXL26, +Sytox_.0001.fcs", truncate_max_range= F)@exprs)
dat <- dat[`SSC-A`>0 & `FSC-A`>0]

par(mfrow=c(2,2))
# gate1
plot(log10(dat$`FSC-A`), 
     log10(dat$`SSC-A`),
     cex= 0.4, 
     pch= 16,
     col= adjustcolor("grey", 0.1),
     xlim= c(1.25, 4),
     ylim= c(1, 4.5),
     xlab= "FSC-A",
     ylab= "SSC-A")

gate <- matrix(c(2.3,2.3,3.75,3.5,1.5,3,4.5,3), ncol= 2)
polygon(gate[,1], gate[,2])
overlaps <- over(SpatialPoints(log10(as.matrix(dat[, .(`FSC-A`, `SSC-A`)]))),
                 SpatialPolygons(list(Polygons(list(Polygon(gate)), 1))))
dat <- dat[overlaps==1 & !is.na(overlaps)]

# gate2
plot(log10(dat$`FSC-W`), 
     log10(dat$`FSC-A`),
     cex= 0.4, 
     pch= 16,
     col= adjustcolor("grey", 0.1),
     xlim= c(2.5,3.5),
     ylim= c(2.2,4),
     xlab= "FSC-W",
     ylab= "FSC-A")

gate <- matrix(c(2.8,2.8,3,3,2.3,3.45,3.45,2.3), ncol= 2)
polygon(gate[,1], gate[,2])
overlaps <- over(SpatialPoints(log10(as.matrix(dat[, .(`FSC-W`, `FSC-A`)]))),
                 SpatialPolygons(list(Polygons(list(Polygon(gate)), 1))))
dat <- dat[overlaps==1 & !is.na(overlaps)]

# gate2
plot(log10(dat$`B1-A`), 
     log10(dat$`FSC-A`),
     cex= 0.4, 
     pch= 16,
     col= adjustcolor("grey", 0.1),
     xlim= c(2.6, 3.1),
     xlab= "GFP",
     ylab= "FSC-A")

gate <- matrix(c(2.8,2.95,3.05,2.95,2.3,3.4,3.4,2.3), ncol= 2)
polygon(gate[,1], gate[,2])
text(3.25, 2.8, "GFP+")




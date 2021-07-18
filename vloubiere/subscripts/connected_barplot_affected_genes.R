setwd("D:/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(ggalluvial)

dat <- data.table(file= list.files("db/FC_tables/", "PH18_vs_RNA_W18|PHD11_vs_RNA_WKD|PHD9_vs_RNA_WKD|PH29_vs_RNA_W29", 
                                   full.names = T))
dat[, cdition:= gsub("RNA_|epiCancer_Parreno_|.txt", "", basename(file)), file]
dat <- dat[, fread(file), (dat)]
dat[padj>=0.01 | is.na(padj), class:= "unaffected"]
dat[padj<0.01 & log2FoldChange>0, class:= "up"]
dat[padj<0.01 & log2FoldChange<0, class:= "down"]
dat[, keep:= ifelse(all(class=="unaffected"), F, T), V1]
dat <- dat[(keep)]

res <- dcast(dat, V1~cdition, value.var = "class")

#-----------------------------------#
# PLOT
#-----------------------------------#
pdf("pdf/alluvial_plot_timecourse.pdf", width = 10)
plot.new()

#### BLOCS --------------------####
blocs <- melt(res[, V1:PHD9_vs_WKD], id.vars = "V1")
blocs[, variable:= factor(variable, 
                          levels = c("PH18_vs_W18", 
                                     "PHD9_vs_WKD",
                                     "PHD11_vs_WKD",
                                     "PH29_vs_W29"))]
setorderv(blocs, c("variable", "value"))
blocs <- blocs[, .N, .(variable, value)]
blocs[, xleft:= (.GRP*2-2)*1/7, variable]
blocs[, xright:= xleft+1/7]
blocs[, ytop:= cumsum(N)/sum(N), variable]
blocs[, ybottom:= c(0, cumsum(N)/sum(N))[-(.N+1)], variable]
blocs[, Cc:= c("cornflowerblue",
               "lightgrey",
               "tomato")[.GRP], value]
blocs[, rect(xleft[1], 
             ybottom[1], 
             xright[1], 
             ytop[1], 
             col= Cc[1]), (blocs)]
blocs[, text(mean(c(xleft, xright)),
             mean(c(ytop, ybottom)), 
             paste0(value, "\n(", N, ")")), (blocs)]
blocs[, text(mean(c(xleft, xright)), 
             1.025,
             variable,
             pos= 3, 
             xpd= T), .(xleft, xright, variable)]


#### Transitions --------------------####
links <- copy(res[, .(PH18_vs_W18, 
                      PHD9_vs_WKD, 
                      PHD11_vs_WKD, 
                      PH29_vs_W29)])
setorderv(links, c("PH18_vs_W18", "PHD9_vs_WKD"))
links[, y1_left:= .I/.N]
setorderv(links, c("PHD9_vs_WKD", "PH18_vs_W18"))
links[, y1_right:= .I/.N]
setorderv(links, c("PHD9_vs_WKD", "PHD11_vs_WKD"))
links[, y2_left:= .I/.N]
setorderv(links, c("PHD11_vs_WKD", "PHD9_vs_WKD"))
links[, y2_right:= .I/.N]
setorderv(links, c("PHD11_vs_WKD", "PH29_vs_W29"))
links[, y3_left:= .I/.N]
setorderv(links, c("PH29_vs_W29", "PHD11_vs_WKD"))
links[, y3_right:= .I/.N]

pl <- rbindlist(list(links[, .(V1= PH18_vs_W18, V2= PHD9_vs_WKD, yleft= y1_left, yright= y1_right)],
                     links[, .(V1= PHD9_vs_WKD, V2= PHD11_vs_WKD, yleft= y2_left, yright= y2_right)],
                     links[, .(V1= PHD11_vs_WKD, V2= PH29_vs_W29, yleft= y3_left, yright= y3_right)]), 
                idcol = T)
pl[, xleft:= (.GRP*2-2)*1/7+1/7, .id]
pl[, xright:= xleft+1/7]

selCc <- function(V1, V2)
{
   if(V1==V2)
      "grey" else if(V1=="up" & V2 %in% c("unaffected", "down"))
         "cornflowerblue" else if(V1=="unaffected" & V2=="down")
            "cornflowerblue" else if(V1=="down" & V2 %in% c("unaffected", "up"))
               "tomato" else if(V1=="unaffected" & V2=="up")
                  "tomato" 
}
pl[, {
   polygon(c(xleft, 
             xright, 
             xright, 
             xleft),
           c(max(yleft), 
             max(yright), 
             min(yright), 
             min(yleft)),
           col= adjustcolor(selCc(V1, V2), 0.5), 
           border= NA)
}, .(.id, V1, V2, xleft, xright)]

dev.off()
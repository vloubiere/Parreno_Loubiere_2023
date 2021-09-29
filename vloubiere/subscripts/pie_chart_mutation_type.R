dat <- readRDS("Rdata/DNA_novogene_clean_final.rds")

for(class in c("INDEL", "SNP"))
{
  pl <- data.table(type= unique(unlist(strsplit(dat[CLASS==class, ANNO_REGION], ","))))
  pl[, Cc:= matlab.like(.N)]
  all <- dat[CLASS==class & variable=="ph18", ANNO_REGION]
  pl[as.data.table(table(unlist(strsplit(all, ",")))), 
     all:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="ph29" & !(loh_from_ph18), ANNO_REGION], ",")))), 
     ph29:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="phd11" & !(loh_from_ph18), ANNO_REGION], ",")))), 
     phd11:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="ph29_t5" & !(loh_from_ph18), ANNO_REGION], ",")))), 
     ph29_t5:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="phd11_t8" & !(loh_from_ph18), ANNO_REGION], ",")))), 
     phd11_t8:= i.N, 
     on= "type==V1"]
  
  cols <- colnames(pl)[-c(1,2)]
  pl[, (cols):= lapply(.SD, function(x) as.numeric(ifelse(is.na(x), 0 , x))), .SDcols= cols]
  pl <- melt(pl, id.vars = c("type", "Cc"))
  
  pdf(paste0("pdf/pie_charts_", class, "s_function.pdf"), 7, 10)
  par(mfrow=c(3,2), 
      mar= c(5,5,5,5))
  pl[, {
    pie(value, 
        col= Cc, 
        labels = type, 
        main= paste0(class, " ", variable, " (", sum(value), ")"),
        cex= 0.8)
    }, by= variable]
  dev.off()
  
  # EXONS
  pl <- data.table(type= unique(unlist(strsplit(dat[CLASS==class, EXON_TYPE_of_MUTATION], ","))))
  pl[, Cc:= matlab.like(.N)]
  all <- dat[CLASS==class & variable=="ph18", EXON_TYPE_of_MUTATION]
  pl[as.data.table(table(unlist(strsplit(all, ",")))), 
     all:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="ph29" & !(loh_from_ph18), EXON_TYPE_of_MUTATION], ",")))), 
     ph29:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="phd11" & !(loh_from_ph18), EXON_TYPE_of_MUTATION], ",")))), 
     phd11:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="ph29_t5" & !(loh_from_ph18), EXON_TYPE_of_MUTATION], ",")))), 
     ph29_t5:= i.N, 
     on= "type==V1"]
  pl[as.data.table(table(unlist(strsplit(dat[CLASS==class & variable=="phd11_t8" & !(loh_from_ph18), EXON_TYPE_of_MUTATION], ",")))), 
     phd11_t8:= i.N, 
     on= "type==V1"]
  
  cols <- colnames(pl)[-c(1,2)]
  pl[, (cols):= lapply(.SD, function(x) as.numeric(ifelse(is.na(x), 0 , x))), .SDcols= cols]
  pl <- melt(pl, id.vars = c("type", "Cc"))
  
  pdf(paste0("pdf/pie_charts_exon_", class, "s_function.pdf"), 7, 10)
  par(mfrow=c(3,2), 
      mar= c(5,5,5,5))
  pl[, {
    pie(value, 
        col= Cc, 
        labels = type, 
        main= paste0(class, " ", variable, " (", sum(value), ")"),
        cex= 0.8)
  }, by= variable]
  dev.off()
}

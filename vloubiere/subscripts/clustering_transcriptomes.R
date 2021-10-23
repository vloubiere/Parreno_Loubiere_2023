setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)

dat <- readRDS("Rdata/final_FC_table.rds")
# Define groups
dat[grepl("dose$", cdition), group:= "dose"]
dat[grepl("TS$", cdition), group:= "epiCancer"]
dat[, group:= factor(group, levels = c("dose", "epiCancer"))]
# Define subgroups
dat[, subgroup:= group]
dat[group=="epiCancer", subgroup:= ifelse(grepl("^EZ", cdition), "EZ", "PH")]
# Define genes to use for clustering
dat[, clustering:= any(padj<=0.05 
                       & abs(log2FoldChange)>1), .(group, FBgn)]
# epiCancer: remove genes affected at 18C 
dat[group=="epiCancer" & clustering, 
    clustering:= all(padj[grepl("18_TS$", cdition)] > 0.05), FBgn]

#---------------------------------------------#
#### Alluvial plots ####
#---------------------------------------------#
pdf("pdf/alluvial_plots.pdf", 
    width = 16/9*7, 
    height = 10)
dat[, {
  mat <- .SD[, {
    FC <- rep("stable", .N)
    FC[padj<=0.05 & log2FoldChange>1] <- "up"
    FC[padj<=0.05 & log2FoldChange<1] <- "down"
    .(FBgn, FC)
  }, cdition]
  mat <- dcast(mat, 
               FBgn~cdition, 
               value.var = "FC")
  mat <- na.omit(mat)[,-1]
  vl_alluvial_plot(mat, 
                   class_levels = c("down", "stable", "up"))
}, subgroup]
dev.off()

#---------------------------------------------#
#### clustering ####
#---------------------------------------------#
pdf("pdf/clustering/heatmaps_clustering.pdf", 
    width = 4)
dat[, cl:= {
  .c <- .SD[, {
    x <- log2FoldChange
    lim <- quantile(x, c(0.005, 0.995))
    x[x<lim[1]] <- lim[1]
    x[x>lim[2]] <- lim[2]
    .(FBgn, log2FoldChange, logFC_clip= x)
  }, cdition]
  # mat for clustering
  mat <- dcast(.c[(clustering)],
               FBgn~cdition,
               value.var = "logFC_clip")
  mat <- as.matrix(mat, 1)
  # mat for plotting
  dmat <- dcast(.c[(clustering)],
                FBgn~cdition,
                value.var = "log2FoldChange")
  dmat <- as.matrix(dmat, 1)
  # Cluster using kmeans
  pl <- vl_heatmap(mat,
                   newdata = dmat,
                   cluster_cols = F,
                   main = group,
                   breaks= c(-3,0,3),
                   clustering_distance_rows = "euclidean",
                   kmeans_k = switch(as.character(group),
                                     "dose"= 6,
                                     "epiCancer"= 9),
                   show_rownames = F, 
                   legend_title = "log2FoldChange")
  rect(0,0,1,1,lwd=0.25, xpd=T)
  pl[, text(0,
            mean(y)/max(pl$y)-grconvertY(1, "chars", "nfc")/2,
            paste0(length(unique(row)), "\n", "genes"),
            pos= 2,
            xpd= T,
            cex= 0.8), rcl]
  pl[.SD, rcl, on= c("row==FBgn", "col==cdition")]
}, group]
dev.off()

#-------------------------------------------------#
#### Heatmaps selected genes ####
#-------------------------------------------------#
genes <- as.data.table(readxl::read_excel("Rdata/list_genes_interest.xlsx"))
genes[, class:= factor(group, 
                       levels= c("PRC1", "PRC2", "HOX", "JAK-STAT", "Cancer markers", 
                                 "RDGN", "Ecdysone metamorhphosis", "cell-cell adhesion GO"))]
genes$group <- NULL

pdf("pdf/hypothesis_driven_heatmaps_epiCancer.pdf", 
    height = 16*0.8, 
    width = 9*0.8)
par(xaxs= "i",
    yaxs= "i")
layout(matrix(c(1,2,3,4,5,6,7,8,8,8,8,8,8,8), 
              ncol= 2), 
       heights = c(6/36, 5/36, 11/36, 4/36, 5/36, 6/36, 9/36), 
       widths = c(0.8,1))
dat[, {
  sub <- dcast(.SD, 
               symbol~cdition, 
               value.var = "log2FoldChange")
  setkeyv(sub, "symbol")
  
  # PLOT
  genes[, {
    mat <- as.matrix(na.omit(sub[symbol]), 1)
    rownames(mat) <- paste0(rownames(mat), "  ")
    colnames(mat) <- paste0(colnames(mat), "  ")
    par(mai= c(ifelse(class %in% c("Ecdysone metamorhphosis",
                                   "cell-cell adhesion GO"), 0.8, 0.1),
               0.75,
               0.3,
               ifelse(class=="cell-cell adhesion GO", 1, 0.15)))
    vl_heatmap(mat,  
               cluster_cols = F,
               cluster_rows = T,
               show_colnames = ifelse(class %in% c("Ecdysone metamorhphosis",
                                                   "cell-cell adhesion GO"), T, F),
               legend_title = "log2FC",
               auto_margins = F, 
               breaks = c(-5,0,5), 
               display_numbers = T)
    box(lwd= 0.25)
    title(class, line = 1)
  }, keyby= class]
}, group]
dev.off()

#---------------------------------------------#
#### Heatmap apsptosis genes ####
#---------------------------------------------#
apop <- dat[grepl("dose", cdition) & 
              symbol %in% c("hid", "rpr", "Diap1", "Drice", "DCP1", "Decay", "Mmp1", "Dronc", "Dredd", "p53")]
apop <- dcast(apop, 
              symbol~cdition, 
              value.var = "log2FoldChange")
apop <- as.matrix(apop, 1)

pdf("pdf/apoptosis_heatmap.pdf", 
    width = 4, 
    height = 7)
vl_heatmap(apop,
           display_numbers = T, 
           legend_title = "log2FoldChange", 
           cluster_cols = F)
dev.off()


#---------------------------------------------#
#### GOs ####
#---------------------------------------------#
pdf("pdf/clustering/GOs_clustering.pdf", 
    width = 9.5, 
    height = 15)
dat[, {
  .l <- split(FBgn, cl)
  .l <- lapply(.l, unique)
  vl_GO_clusters(FBgn_list = .l,
                 padj_cutoff =  0.0001,
                 N_top = switch(as.character(group),
                                "dose"= Inf,
                                "epiCancer"= 11),
                 all_FBgns = unique(FBgn),
                 cex = 0.3)
  print("")
}, group]
dev.off()

#---------------------------------------------#
#### STRING ####
#---------------------------------------------#
cols <- c("pink","pink3","indianred2","red2",
          "darkorchid1","purple","darkorchid4",
          "olivedrab1","limegreen","olivedrab3","olivedrab4",
          "lightsteelblue1","cornflowerblue","blue","navy",
          "black","gold","goldenrod","goldenrod4","sienna4")
dat[, Cc:= {
  Cc <- colorRampPalette(cols)(max(cl, na.rm=T))
  Cc[cl]
}, group]

pdf("pdf/clustering/STRING_networks_clustering.pdf", 
    width = 20, 
    height = 20)
dat[, {
  .c <- .SD[, .(value= max(abs(log2FoldChange))), .(symbol, cl, Cc)]
  .c <- na.omit(unique(.c))
  clip_size <- .c$value*0.7
  clip_size[clip_size>10] <- 10
  clip_size[clip_size<2] <- 2
  obj <- vl_STRING_interaction(symbols = .c$symbol,
                               size = clip_size, 
                               cex.label = clip_size,
                               col = adjustcolor(.c$Cc, 0.6),
                               score_cutoff = switch(as.character(group),
                                                     "epiCancer" = 800,
                                                     "dose"= 900))
  set.seed(123)
  vl_STRING_network(obj, 
                    cex.vertices.labels = 0.2, 
                    vertex.border.col = NA)
  leg <- unique(.c[, .(cl, Cc)][order(cl)])
  legend("topleft", 
         bty= "n",
         legend = paste0("cluster ", leg$cl),
         fill= leg$Cc,
         cex= 2)
  print("")
}, group]
dev.off()

#---------------------------------------------#
#### Motifs ####
#---------------------------------------------#
pdf("pdf/clustering/motif_enrichment_clustering.pdf", 
    width = 9, 
    height = 9)
par(mfrow=c(1,2))
dat[, {  
  file <- paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/motif_enrichment_tables/", 
                 group, "_motif_enrichment_between_clusters.rds") #Using genes form other clusters as background
  if(!file.exists(file))
  {
    #### Import motif counts table and add cl column
    if(!exists("mot"))
      mot <- readRDS("Rdata/ED_REs_motif_counts_final.rds")
    if("cl" %in% names(mot))
      mot$cl <- NULL
    mot[unique(.SD[, .(FBgn, cl)]), cl:= i.cl, on= "FBgn"]
    #### Compute enrichment
    obj <- vl_motif_cl_enrich(obj = mot, 
                              bg = na.omit(unique(mot$cl)), 
                              cl_column = "cl")
    saveRDS(obj, file)
  }else
  {
    obj <- readRDS(file)
    # Select top enrich for each TF
    obj <- obj[motif %in% obj[, motif[which.min(padj)], motif_name]$V1]
  }
  
  # PLOT
  par(mar= c(5, switch(as.character(group),
                       "dose"= 12,
                       "epiCancer"= 10), 2, 5), 
      cex=1)
  pl <- vl_motif_cl_enrich_plot_only(obj, 
                                     padj_cutoff = 0.001, 
                                     N_top = switch(as.character(group),
                                                    "dose"= 10,
                                                    "epiCancer"= 6), 
                                     auto_margin = F)
  # heatmap FC TFs 
  mat <- dcast(.SD,
               symbol~cdition,
               value.var = "log2FoldChange")
  mat <- as.matrix(mat, 1)
  sel <- unique(pl[!is.na(cl)][order(y, decreasing = T), motif_name])
  sel <- unlist(tstrsplit(sel, "/", keep= 1))
  mat <- mat[match(sel, rownames(mat)),]
  par(mar= c(5, switch(as.character(group),
                       "dose"= 14,
                       "epiCancer"= 14), 2, 5),
      cex= 0.7)
  pl <- vl_heatmap(mat,
                   cluster_rows = F,
                   cluster_cols = F,
                   breaks = c(-3, 0, 3),
                   auto_margin = F, 
                   show_rownames = F, 
                   legend_title = "log2FC", 
                   display_numbers = T)
  rect(0, 0, 1, 1, lwd=0.25)
  pl[, text(0, 
            y/max(pl$y)-(1/max(pl$y)/2), 
            unlist(tstrsplit(row, "_", keep= 1)), 
            pos= 2, 
            xpd= T), .(y, row)]
  print("")
}, group]
dev.off()

#---------------------------------------------#
#### Make excel file ####
#---------------------------------------------#
df <- dcast(dat, 
            FBgn+symbol~cdition, 
            value.var = list("log2FoldChange", "padj"))
# Add cluster
dat[, {
  col <- paste0(group, "_cl")
  df[unique(.SD[, .(FBgn, cl)]), (col):= i.cl, on="FBgn"]
  print("")
}, group]
setcolorder(df, c("FBgn", "symbol", "dose_cl", "epiCancer_cl"))

# Rename columns
names(df)[grep("^log2FoldChange", names(df))] <- gsub("^log2FoldChange_(.*)", "\\1_log2FC", names(df)[grep("^log2FoldChange", names(df))])
names(df)[grep("^padj", names(df))] <- gsub("^padj_(.*)", "\\1_padj", names(df)[grep("^padj", names(df))])

# create a workbook
require(openxlsx)
wb <- openxlsx::createWorkbook() 
#add a worksheet to the workbook
openxlsx::addWorksheet(wb, "Sheet", 
                       gridLines = TRUE) 
# write my analysis into the worksheet of the workbook
openxlsx::writeData(wb, 
                    "Sheet", 
                    df, 
                    keepNA = T) 
# Create a list of cell format depending on padj/FC
Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))(7)
checks <- list(down3= list(test= "AND($col11<0,$col21<0.05)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[3])),
               down2= list(test= "AND($col11<0,$col21<0.01)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[2])),
               down1= list(test= "AND($col11<0,$col21<0.001)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[1])),
               down1= list(test= "AND($col11>0,$col21<0.05)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[5])),
               down2= list(test= "AND($col11>0,$col21<0.01)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[6])),
               down3= list(test= "AND($col11>0,$col21<0.001)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[7])))
# Create canonical excel col names
columns <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
columns <- c(columns, paste0("A", columns))
columns <- c(columns, paste0("B", columns))
# account padj and FC to colour cells
for(icol in grep("_log2FC$", names(df)))
{
  for(j in seq(checks))
  {
    conditionalFormatting(wb, 
                          "Sheet", 
                          cols = icol,
                          rows = 1:nrow(df),
                          rule = gsub("col2", columns[icol+12], gsub("col1", columns[icol], checks[[j]]$test)), 
                          style = checks[[j]]$style,
                          type = "expression")
  }
}
dir.create("tables_AMM", 
           showWarnings = F)
saveWorkbook(wb, 
             "tables_AMM/final_table.xlsx", 
             overwrite = TRUE)

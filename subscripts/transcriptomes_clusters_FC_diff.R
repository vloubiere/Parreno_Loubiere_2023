require(data.table)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[is.na(cl), cl:= 0]
# Select columns to plot
cols <- c("H3K27me3_PH18_body",
          "H3K27me3_PH29_body",
          "H3K27me3_PHD11_body",
          "H2AK118Ub_PH18_body",
          "H2AK118Ub_PH29_body",
          "H2AK118Ub_PHD11_body")
# Melt data
.m <- melt(dat[, c("cl", "PRC1_bound", cols), with= F], 
           id.vars = c("cl", "PRC1_bound"), 
           measure.vars = cols)
.m[, ylim:= fcase(grepl("K118Ub", variable), 15,
                  grepl("K27me3", variable), 30)]
.m[, cdition:= fcase(grepl("_PH18_", variable), "No ph-KD",
                     grepl("_PH29_", variable), "Constant ph-KD",
                     grepl("_PHD11_", variable), "Transient ph-KD")]

##############################################
# PLOT
##############################################
pdf("pdf/cluster_PRC1_bound_unbound_RNA.pdf",
    width= 6,
    height= 3.5)
par(las= 2,
    mar= c(3.5,3.5,2,0.5),
    mgp= c(2, 0.5, 0),
    tcl= -0.2,
    mfrow= c(2,3),
    lwd= 0.5)
.m[, {
  box <- vl_boxplot(value~PRC1_bound+cl,
                    .SD,
                    col= c("lightgrey", "tomato"),
                    compute_pval = list(c(1,2), c(3,4), c(5,6), c(7,8), c(9,10), c(11,12), c(13,14)),
                    xaxt= "n",
                    at= rep(seq(1, 13, 2), each= 2)+c(0.2, 0.8),
                    ylim= c(0, ylim),
                    ylab= paste(tstrsplit(variable, "_", keep= 1), "enrichment"),
                    main= cdition)
  text(x= rep(box$x, each= 2)+strwidth("M", cex= 0.4)*c(-1,1), 
       y= rep(box$y, each= 2)+strheight("M", cex= 0.5),
       labels= paste0("n= ", box$n),
       col= c("lightgrey", "tomato"),
       pos= 4,
       srt= 90,
       offset= 0,
       cex= 0.5)
  axis(1, at= seq(1.5, 13.5, 2), labels = NA)
  vl_tilt_xaxis(x = seq(1.5, 13.5, 2),
                labels =  c("Unaffected", paste("Cluster", 1:6)),
                srt= 30)
  legend("topleft",
         legend= c("PRC1+", "PRC1-"),
         fill= c("tomato", "lightgrey"),
         bty= "n",
         xpd=T,
         cex= 0.8)
  abline(h= box$stats[3,1], lty= 2)
  print("done")
}, .(variable, cdition, ylim)]
dev.off()

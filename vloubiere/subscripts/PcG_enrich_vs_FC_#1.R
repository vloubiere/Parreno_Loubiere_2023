setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

cl <- readRDS("Rdata/clustering_RNA.rds")
obj <- cl$data[abs(RNA_PH29)>1 & RNA_PH29_padj<0.05]
tracks <- c("db/bw/SA_2020/PC_ED_merge.bw",
            "db/bw/SA_2020/PH_ED_merge.bw",
            "db/bw/SA_2020/PSC_ED_merge.bw",
            "db/bw/SA_2020/H2AK118Ub_ED_merge.bw",
            "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
            "db/bw/cutnrun/H3K27Ac_PH29_merge.bw",
            "db/bw/cutnrun/H3K27Ac_PHD11_merge.bw",
            "db/bw/cutnrun/H3K27Ac_PHD9_merge.bw",
            "db/bw/cutnrun/H3K27me3_PH18_merge.bw",
            "db/bw/cutnrun/H3K27me3_PH29_merge.bw",
            "db/bw/cutnrun/H3K27me3_PHD11_merge.bw",
            "db/bw/cutnrun/H3K27me3_PHD9_merge.bw")

pdf("pdf/Figures/PcG_enrich_revert_vs_not_#1.pdf",
    height = 17.5, 
    width = 60)
layout(matrix(1:(6*24), nrow=6, byrow = T), 
       widths= rep(c(1,0.5), 12))
for(cdition in c("PHD9", "PHD11")) # Use the two conditions to split recovering genes
{
  for(cluster in list(c(2,5), 2,5)) # Each clsuter separately or merge
  {
    for(track in tracks) # Each bw track
    {
      dat <- obj[cl %in% cluster]
      # Split genes into classes
      if("class" %in% names(dat))
        dat$class <- NULL
      if(cdition=="PHD9")
      {
        dat[RNA_PHD9>1 & RNA_PHD9_padj<0.05, class:= "UP D9"]
        dat[is.na(class), class:= "NOT UP D9"]
      }
      if(cdition=="PHD11")
      {
        dat[RNA_PHD11>1 & RNA_PHD11_padj<0.05, class:= "UP D11"]
        dat[is.na(class), class:= "NOT UP D11"]
      }
      dat[, class:= paste0(ifelse(is.na(PRC1_cluster), "PRC1-", "PRC1+"), ": UP 29 , ", class)]
      dat[, class:= factor(class, rev(unique(sort(class))))]
      # Plot average track for each class
      par(mar= c(5,4,2,1),
          las= 1,
          mgp= c(2.5,0.5,0),
          tcl= -0.2)
      .q <- vl_bw_average_track(bed= GRanges(dat[, TSS]), 
                                tracks= track,
                                set_IDs = dat[, class],
                                upstream = 2500,
                                downstream = 10000,
                                stranded = T,
                                center_label = "TSS", 
                                legend.cex = 0.6, 
                                legend= F)
      # I bw track is for a protein, keep only TSS-prox bins
      if(grepl("PH_ED|PC_ED|PSC_ED", track))
      {
        .q <- .q[between(bin.x, 13, 26, incbounds = T)]
        abline(v= c(13, 26), lty=2)
        xl <- "[quantification window]"
        text(26, 
             grconvertY(0.75, "npc", "user")+par("usr")[3],
             "Quantification window",
             pos= 4,
             cex= 0.6)
      } else
        xl <- "[whole window]"
      # Titlte label
      title(paste0(gsub("_merge.bw", "", basename(track)), 
                   " cluster ", 
                   paste0(cluster, collapse = "+"),
                   "  / ", cdition))
      # Quantification boxplot
      par(mar= c(5,4,6,2))
      x <- .q[, .(score= mean(score)), .(set_IDs, region_ID)]
      x <- split(x$score, x[, set_IDs], drop= T)
      vl_boxplot(x,
                 violcol= .q[, col[1], keyby= set_IDs]$V1,
                 violin=T, 
                 xaxt= "n",
                 ylab= "Enrichment", 
                 compute_pval= list(c(1,2)))
      # Legend
      title(xlab= xl)
      leg <- .q[, .(col= col[1]), keyby= set_IDs]
      leg[dat, N:= .N, set_IDs, on= "set_IDs==class"]
      legend(grconvertX(0, "nfc", "user"),
             grconvertY(1, "nfc", "user")-diff(grconvertY(c(0,2), "lines", "user")),
             xpd= T, 
             legend= leg[, paste0(set_IDs, " (", N, ")")],
             fill= leg$col,
             bty= "n",
             cex= .7)
    }
  }
}
dev.off()

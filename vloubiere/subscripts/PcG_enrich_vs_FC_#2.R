setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(GenomicRanges)
require(vlfunctions)

# Import genes which are up-regulated at 29
dat <- fread("db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PH29_vs_RNA_W29.txt")
dat <- dat[padj<0.05 & log2FoldChange>1]
# Add extended promoters
proms <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
proms <- as.data.table(proms)
dat[proms[type=="gene"], c("symbol", "seqnames", "start", "end", "strand"):= .(i.gene_symbol, paste0("chr", i.seqnames), i.start, i.end, i.strand), on="FBgn==gene_id"]
dat <- vl_resizeBed(dat, "start", 2500, 2500, ignore.strand = F)
# Add PRC1 binding
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
PRC1 <- loadRData("external_data/SA2020_cl.list")
PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
dat[PRC1, PRC1_cluster:= i.PRC1_cluster, on= "symbol"]
dat <- dat[!is.na(PRC1_cluster)]
# Overlap K27me3 and PRC1 binding
dat <- dat[vl_covBed(dat, "db/peaks/K27_cutnrun/H3K27me3_PH18_confident_peaks.bed")>0]

# tracks
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

pdf("pdf/Figures/PcG_enrich_revert_vs_not_using_PH29_as_control.pdf",
    height = 17.5, 
    width = 60)
layout(matrix(1:(6*24), nrow=6, byrow = T), 
       widths= rep(c(1,0.5), 12))
for(cdition in c("PHD9", "PHD11"))
{
  # Import PHDX vs PH29 foldChanges
  if(cdition=="PHD9")
    recov <- fread("db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD9_vs_RNA_PH29.txt") else if(cdition=="PHD11")
      recov <- fread("db/FC_tables/RNA/epiCancer_ED_GFP-_system_RNA_RNA_PHD11_vs_RNA_PH29.txt")
  # Draw classes
  if("class" %in% names(dat))
    dat$class <- NULL
  dat[FBgn %in% recov[padj<0.05 & log2FoldChange<(-1), FBgn], class:= "PRC1+: UP 29 , RECOVERY"]
  dat[is.na(class), class:= "PRC1+: UP 29 , NO RECOVERY"]
  for(track in tracks) # Each bw track
  {
    # Plot average track for each class
    par(mar= c(5,4,2,1),
        las= 1,
        mgp= c(2.5,0.5,0),
        tcl= -0.2)
    .q <- vl_bw_average_track(bed= dat, 
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
    title(paste0(gsub("_merge.bw", "", basename(track)), " Using ", cdition, "_vs_PH29"))
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
dev.off()

setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

# Import
dat <- fread("Rdata/RECOVERY_NORECOVERY_genes.txt")
dat <- dat[!is.na(PHD9_RECOVERY) & !is.na(PHD11_RECOVERY) & PHD9_RECOVERY==PHD11_RECOVERY]
dat <- dat[, .(FBgn, symbol, recov= PHD11_RECOVERY)]
# Add gene coor
genes <- rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[type=="gene", .(FBgn= gene_id, seqnames, start, end)]
dat <- genes[dat, on= "FBgn"]
# Tracks
bw <- c("db/bw/cutnrun/H3K27me3_PH18_merge.bw",
        "db/bw/cutnrun/H3K27Ac_PH18_merge.bw",
        "db/bw/cutnrun/H3K4me1_PH18_merge.bw",
        "db/bw/RNA-Seq/RNA_PH18_merge.bw")
bw2 <- c("db/bw/cutnrun/H3K27me3_PHD11_merge.bw",
         "db/bw/cutnrun/H3K27Ac_PHD11_merge.bw",
         "db/bw/cutnrun/H3K4me1_PHD11_merge.bw",
         "db/bw/RNA-Seq/RNA_PHD11_merge.bw")
bw3 <- c("db/bw/cutnrun/H3K27me3_PH29_merge.bw",
         "db/bw/cutnrun/H3K27Ac_PH29_merge.bw",
         "db/bw/cutnrun/H3K4me1_PH29_merge.bw",
         "db/bw/RNA-Seq/RNA_PH29_merge.bw")
cols1 <- adjustcolor("blue", 0.25)
cols2 <- adjustcolor("green", 0.25)
cols3 <- adjustcolor("red", 0.25)
max <- c(30,30,30,50)
min <- c(0,0,0,5)

pdf("pdf/Figures/RECOVERY_vs_NOT_screenshot.pdf", width = 10)
dat[c("Abd-B", "Psc", "upd1", "dpp"), {
  vl_screenshot(bed= data.table(seqnames, start= start-10000, end= end+10000),
                tracks= bw,
                names= c("H3K27me3", "H3K27Ac", "H3K4me1", "RNA"),
                max= max,
                min= min,
                col= cols1,
                genome = "dm6")
  vl_screenshot(bed= data.table(seqnames, start= start-10000, end= end+10000),
                tracks= bw2,
                max= max,
                min= min,
                add= T, 
                bg_col = adjustcolor("white", 0),
                col= cols2)
  vl_screenshot(bed= data.table(seqnames, start= start-10000, end= end+10000),
                tracks= bw3,
                max= max,
                min= min,
                add= T, 
                bg_col = adjustcolor("white", 0),
                col= cols3)
  title(main= ifelse(recov, "RECOVERY", "NO RECOVERY"))
}, on= "symbol"]
dev.off()
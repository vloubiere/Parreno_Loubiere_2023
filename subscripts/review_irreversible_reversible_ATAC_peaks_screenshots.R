setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)

# Regions to plot ----
sel <- data.table(ID= c("chrX:18249640-18299430",
                        "chr3R:16656000-16710406"))
                        # "chr3R:16741987-16758726"))
sel[, c("seqnames", "start", "end"):= vl_toDTranges(ID)]

# Save bed files ----
dat <- readRDS("Rdata/final_ATAC_table.rds")
file.remove(list.files("db/bed/ATAC_peaks/cluster/", full.names = T))
meta <- dat[!is.na(cl), {
  file <- paste0("db/bed/ATAC_peaks/cluster/", cl, "_peaks.bed")
  fwrite(.SD[, .(seqnames, start, end)],
         file,
         col.names = F,
         row.names = F,
         sep= "\t",
         quote= F)
  .(file= file)
}, keyby= cl]
meta <- meta[-1]
meta[, cl:= droplevels(cl)]

# bw tracks ---
tracks <- c("db/bw/ATAC/ATAC_PH18_merge.bw",
            "db/bw/ATAC/ATAC_PH29_merge.bw",
            "db/bw/ATAC/ATAC_PHD11_merge.bw",
            meta$file)
names <- c("No ph-KD", "Constant ph-KD", "Transient ph-KD",
           levels(meta$cl))

# Plot ----
pdf("pdf/review_screenshot_ATAC_reversible_irreversible_genes.pdf",
    width = 5, 
    height = 2.5)
vl_par(mai= c(.7, 1.2, .3, .2),
       mgp= c(0,0,0))
vl_screenshot(bed = sel,
              tracks= tracks, 
              names = names, 
              widths = c(100L, 60L),
              max= c(25,25,25),
              genome = "dm6")
dev.off() 
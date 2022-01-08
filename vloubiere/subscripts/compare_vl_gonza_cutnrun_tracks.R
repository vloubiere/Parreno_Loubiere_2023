setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(GenomicRanges)
require(BSgenome.Dmelanogaster.UCSC.dm6)
require(kohonen)

#-----------------------------------------------------------#
# Screenshots replicates
#-----------------------------------------------------------#
tracks <- data.table(file= list.files(c("db/bw/cutnrun_reps_gonza",
                                        "db/bw/cutnrun_reps_vl"), full.names = T))
tracks[, cdition:= tstrsplit(basename(file), "_rep", keep= 1)]
tracks[, user:= ifelse(grepl("_vl", file), "vl", "gonza")]
tracks <- tracks[, .(seqnames= "chr3R",
                     start= c(16.615e6, 6.625e6, 5e6),
                     end= c(17.025e6, 7.1e6, 15e6)), (tracks)]

pdf("pdf/cutnrun/screenshot_replicates_gonza_vl.pdf", 
    width = 10, 
    height = 5)
tracks[, {
  .n <- paste0(gsub(".bw$", "", basename(file)), "_", user)
  vl_screenshot(GRanges(seqnames,
                        IRanges(start,
                                end)),
                file, 
                names = .n,
                genome = "dm6")
},  .(cdition, seqnames, start, end)]
dev.off()

#-----------------------------------------------------------#
# Screenshots merge
#-----------------------------------------------------------#
tracks <- data.table(file= list.files(c("db/bw/cutnrun_merge_gonza/",
                                        "db/bw/cutnrun_merge_vl/"), full.names = T))
tracks[, cdition:= tstrsplit(basename(file), "_merge", keep= 1)]
tracks[, user:= ifelse(grepl("_vl", file), "vl", "gonza")]
tracks <- tracks[, .(seqnames= "chr3R",
                     start= c(16.615e6, 6.625e6, 5e6),
                     end= c(17.025e6, 7.1e6, 15e6)), (tracks)]

pdf("pdf/cutnrun/screenshot_merge_gonza_vl.pdf",
    width = 10,
    height = 6)
tracks[, {
  .n <- paste0(gsub(".bw$", "", basename(file)), "_", user)
  vl_screenshot(GRanges(seqnames,
                        IRanges(start,
                                end)),
                file, 
                names = .n,
                genome = "dm6")
},  .(cdition, seqnames, start, end)]
dev.off()

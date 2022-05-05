setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)
require(rtracklayer)

# Import
dat <- data.table(file= list.files("db/FC_tables/cutnrun/", full.names = T))
dat[, c("ChIP", "cdition"):= tstrsplit(basename(file), "_", keep= 1:2)]
dat <- dat[, fread(file)[abs(log2FoldChange)>1 & padj<0.05], (dat)]
dat[, state:= ifelse(log2FoldChange>0, "Up", "Down")]

# Compute genes overlaping affected regions
genes <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(genes) <- "UCSC"
genes <- as.data.table(genes)
genes <- genes[type=="gene" & gene_id %in% genes[type %in% c("mRNA", "ncRNA"), gene_id]]
genes <- vl_resizeBed(genes, center = "start", upstream = 2000, downstream = 2000)
ov <- genes[dat, .(FBgn= .(gene_id), symbol= .(gene_symbol)), .EACHI, on= c("seqnames", "start<=end", "end>=start")]
dat[, c("FBgn", "symbol"):= ov[, .(FBgn, symbol)]]
fwrite(dat,
       sep= "\t",
       "db/FC_tables/CUTNRUN_table_AMM.txt")

pdf("pdf/cutnrun/GO_genes_overlapping_CUTNRUN_changes.pdf", 
    width = 9,
    height = 12)
par(las= 2)
dat[, {
  .c <- .SD[, .(FBgn= na.omit(unique(unlist(FBgn)))), keyby= .(cdition= paste0(state, "_", cdition))]
  vl_GO_clusters(split(.c$FBgn, .c$cdition),
                 padj_cutoff= 1e-5, 
                 cex.balloons= 0.5)
}, ChIP]
dev.off()

# Overlap with corresponding transcriptome
ov <- dat[, {
  FC <- fread(list.files("db/FC_tables/RNA/", 
                         paste0("epiCancer_ED_RNA_CUTNRUN_RNA_", cdition, "_vs_RNA_W"),
                         full.names= T))
  FC[, state:= fcase(log2FoldChange>1 & padj<0.05, "Up",
                     log2FoldChange<-1 & padj<0.05, "Down",
                     default= "Unaffected")]
  .c <- .SD[, .(FBgn, state)]
  .c <- na.omit(.c[, .(FBgn= unlist(FBgn)), state])
  .c[FC, transcriptome:= i.state, on= "FBgn"]
  .c[is.na(transcriptome), transcriptome:= "Unaffected"]
}, .(ChIP, cdition)]

pdf("pdf/cutnrun/FBgn_overlaps_cutnrun_transcriptome_changes.pdf",
    width= 15, 
    height = 6)
par(mfrow=c(2,3))
ov[, {
  vl_upset_plot(split(FBgn, 
                      .SD[, .(paste0("ChIP_", state), paste0("RNA_", transcriptome))]))
  title(main= paste(ChIP, cdition))
}, .(ChIP, cdition)]
dev.off()
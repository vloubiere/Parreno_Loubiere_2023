setwd("D:/_R_data/projects/epigenetic_cancer/")
require(Rsubread)
require(DESeq2)
require(DEXSeq)
require(data.table)
require(rtracklayer)
# test
#----------------------------------------------------------#
# Build dm6 index
#----------------------------------------------------------#
# setwd("/_R_data/genomes/dm6/subreadr_index/")
# ref <- "/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
# buildindex(basename= "subreadr_dm6_index", reference= ref)
# setwd("/_R_data/projects/epigenetic_cancer/")

#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
aln <- data.table(file= list.files("db/fastq", ".fq.gz", full.names = T, recursive = T))
aln[, project:= tstrsplit(file, "/", keep= 3)]
aln[, {
  .c <- paste0("db/bam/", project)
  if(!dir.exists(.c)){
    dir.create(.c)
  }}, project]
aln[, bam:= paste0("db/bam/", project, "/", .BY[1], ".bam"), gsub("_1.fq.gz|_2.fq.gz|.fq.gz", "", basename(file))]
aln[, {
  if(!file.exists(bam)){
    print("START...")
    if(.N==2){
      stats <- capture.output(subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", readfile1= file[1], readfile2= file[2], 
                                      maxMismatches = 6, nthreads = 10, unique = T, output_file= bam))
    }else if(.N==1){
      stats <- capture.output(subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", readfile1= file[1], 
                                      nthreads = 10, unique = T, output_file= bam))
    }
    writeLines(stats, con = gsub(".bam$", "_stats.txt", bam))
  }
  print(paste(bam, "DONE!"))
}, bam]

#----------------------------------------------------------#
# Compute counts
#----------------------------------------------------------#
counts <- data.table(file= list.files("db/bam/", ".bam$", recursive = T, full.names = T))
counts[, project:= tstrsplit(file, "/", keep= 4)]
counts[, {
  .c <- paste0("db/counts/", project)
  if(!dir.exists(.c)){
    dir.create(.c)
  }}, project]
counts[, counts_file:= paste0("db/counts/", project, "/", gsub(".bam", "_counts.rds", basename(file)))]
counts[, {
  if(!file.exists(counts_file)){
    if(.N==2){
      counts <- featureCounts(file, annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", isGTFAnnotationFile = T, isPairedEnd = T, nthreads = 8)
    }else if(.N==1){
      counts <- featureCounts(file, annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf", isGTFAnnotationFile = T, isPairedEnd = F, nthreads = 8)
    }
    saveRDS(counts, counts_file)
  }
  print(paste(counts_file, "DONE!"))
}, .(file, counts_file)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dds <- data.table(file= list.files("db/counts/", "_counts.rds$", recursive = T, full.names = T))
dds[, project:= tstrsplit(file, "/", keep= 4)]
dds[, dds_file:= paste0("db/dds/", project, "_dds.rds")]
dds[, rep:= as.numeric(sapply(gsub("_counts.rds", "", .SD$file), function(x) substr(x, nchar(x), nchar(x))))]
.expr <- "_1_counts.rds|_2_counts.rds|_3_counts.rds|_4_counts.rds|_rep1_counts.rds|_rep2_counts.rds|_rep3_counts.rds|_rep4_counts.rds"
dds[, cdition:= gsub(.expr, "", basename(file))]
dds[, cdition:= gsub("_$", "", cdition)]
dds[, {
  sampleTable <- data.frame(unique(.SD[, .(cdition, rep, counts_file= basename(file))]), row.names = "counts_file")
  sampleTable <- apply(sampleTable, 2, as.factor)
  if(!file.exists(dds_file)){
    # DESEq2
    DF <- .SD[, {
      .c <- data.table(readRDS(file)$counts, keep.rownames = T)
      colnames(.c)[2] <- "counts"
      .c
    }, file]
    DF <- data.frame(dcast(DF, rn~basename(file), value.var = "counts"), row.names = "rn")
    DF <- DF[rowSums(DF)>10, ]
    dds <- DESeqDataSetFromMatrix(countData= DF, colData= sampleTable, design= ~rep+cdition)
    dds <- DESeq(dds)
    saveRDS(dds, dds_file)
  }else{
    dds <- readRDS(dds_file)
  }
  print(paste(dds_file, "DONE!"))
  # Differential expression
  diff <- CJ(sampleTable[, "cdition"], sampleTable[, "cdition"], unique = T)
  diff <- diff[V1!=V2]
  diff[, FC_file:= paste0("db/FC_tables/", project, "_", V1, "_vs_", V2, "_FC.txt"), .(V1, V2)]
  diff[, {
    if(!file.exists(FC_file)){
      .c <- as.data.frame(lfcShrink(dds, type= "ashr", contrast= c("cdition", V1, V2)))
      # .c <- as.data.frame(results(dds, contrast= c("cdition", V1, V2)))
      fwrite(.c, FC_file, col.names = T, row.names = T, sep ="\t", quote= F)
    }
    print(paste(FC_file, "DONE!"))
  }, .(FC_file, V1, V2)]
}, .(project, dds_file)]

#----------------------------------------------------------#
# Compute exon counts
#----------------------------------------------------------#
# Create flat exon gtf
# gtf <- as.data.table(import( "../../genomes/dm6/dmel-all-r6.36.gtf"))
# exon <- unique(gtf[type== "exon"][, .(seqnames, start, end, strand, gene_id)])
# exon[, exon_id:= paste0(gene_id, ":", formatC(.SD[, .I], width = 3, flag = "0")), gene_id]
# saf <- data.table(exon[, .(GeneID= exon_id, Chr= paste0("chr", seqnames), Start= start, End= end, Strand= strand)])
# fwrite(saf, "../../genomes/dm6/dm6_r6.36_exons.saf", col.names = T, row.names = F, sep= "\t", quote= F)

if(!file.exists("db/exon_counts/RNA_epiCancer/exon_counts_final_table.txt")){
  files <- list.files("db/bam/RNA_epiCancer/", ".bam$", recursive = T, full.names = T)
  counts <- Rsubread::featureCounts(files, annot.ext= "../../genomes/dm6/dm6_r6.36_exons.saf", isGTFAnnotationFile = F, isPairedEnd = T, nthreads = 8)
  res <- as.data.table(counts$counts, keep.rownames = "ID")
  fwrite(res, "db/exon_counts/RNA_epiCancer/exon_counts_final_table.txt", col.names = T, row.names = F, sep= "\t", quote= F)
}

# make dxd object
if(!file.exists("db/dxd/RNA_epiCancer_dxd.rds")){
  mat <- as.matrix(fread("db/exon_counts/RNA_epiCancer/exon_counts_final_table.txt"), 1)
  mat <- mat[rowSums(mat)>10,]
  
  sampleTable <- data.frame(row.names = colnames(res)[-1], 
                            cdition = unlist(tstrsplit(colnames(res)[-1], "_", keep= 1)),
                            rep = gsub(".bam", "", unlist(tstrsplit(colnames(res)[-1], "_", keep= 2))))
  
  dxd <- DEXSeqDataSet(mat, sampleData = sampleTable, 
                       featureID = as.character(unlist(tstrsplit(res$ID, ":", keep= 2))), 
                       groupID = as.character(unlist(tstrsplit(res$ID, ":", keep= 1))), 
                       design= ~ rep + cdition + exon + cdition:exon,)
  dxd = estimateSizeFactors(dxd)
  dxd = estimateDispersions(dxd)
  saveRDS(dxd, "db/dxd/RNA_epiCancer_dxd.rds")
}
plotDispEsts(dxd)

dxd = testForDEU(dxd, reducedModel = ~ cdition + exon)
dxd = estimateExonFoldChanges(dxd, fitExpToVar= "cdition", denominator = "WKD")
dxr1 = DEXSeqResults( dxd )

mcols(dxr1)$description
table ( dxr1$padj < 0.1 )
table ( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )
plotMA( dxr1, cex=0.8 )

plotDEXSeq(dxr1, "FBgn0000015", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

wh = (dxr2$groupID=="FBgn0010909")
stopifnot(sum(dxr2$padj[wh] < formals(plotDEXSeq)$FDR)==1)

plotDEXSeq( dxr2, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, norCounts=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, splicing=TRUE,
            legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

DEXSeqHTML( dxr2, FDR=0.1, color=c("#FF000080", "#0000FF80") )




















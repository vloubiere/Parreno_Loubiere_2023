setwd("D:/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(Rsubread)
require(DESeq2)

meta <- fread("Rdata/raw_metadata_final.txt")[type=="RNA"]

#----------------------------------------------------------#
# Build dm6 index
#----------------------------------------------------------#
if(!file.exists("D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index.log"))
{
  ref <- "D:/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
  buildindex(basename= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", reference= ref)
}

#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
dir.create("db/bam/", showWarnings = F)
meta[, {
  bam <- paste0("db/bam/", project, "/", prefix_output, ".bam")
  if(!file.exists(bam))
  {
    print("START...")
    subjunc(index= "D:/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index",
            readfile1= paste0("db/fastq/", fq_file[1]),
            readfile2= if(layout=="PAIRED") paste0("db/fastq/", fq_file[2]) else NULL,
            maxMismatches = 6,
            nthreads = 10,
            unique = T,
            output_file= bam)
    print(paste(bam, "--> DONE!"))
  }else
    print(paste(bam, "--> ALREADY EXISTS!"))
}, .(cdition, project, prefix_output, layout)]

#----------------------------------------------------------#
# Counts
#----------------------------------------------------------#
dir.create("db/counts/", showWarnings = F)
meta[, {
  counts <- paste0("db/counts/", project, "/", prefix_output, "_counts.rds")
  if(!file.exists(counts))
  {
    bam <- paste0("db/bam/", project, "/", prefix_output, ".bam")
    .c <- featureCounts(bam,
                        annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf",
                        isGTFAnnotationFile = T,
                        isPairedEnd = ifelse(layout=="PAIRED", T, F),
                        nthreads = 8)
    saveRDS(.c, counts)
    print(paste(counts, "--> DONE!"))
    print(paste(counts, "--> DONE!"))
  }else
    print(paste(counts, "--> ALREADY EXISTS!"))
}, .(project, prefix_output, layout)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dir.create("db/dds/", showWarnings = F)
dds <- meta[type=="RNA"]
dds <- dds[, .(DESeq2_object= unlist(tstrsplit(DESeq2_object, ";"))), setdiff(colnames(dds), "DESeq2_object")]
dds[, 
     {
       dds_file <- paste0("db/dds/", DESeq2_object, "_dds.rds")
       if(!file.exists(dds_file))
       {
         .c <- .SD[, .(file= paste0("db/counts/", project, "/", prefix_output, "_counts.rds")), .(cdition, rep, prefix_output)]
         .c <- unique(.c)
         sampleTable <- data.frame(.c[, .(cdition, rep, prefix_output)], row.names = "prefix_output")
         .c <- .c[, {
           x <- as.data.table(readRDS(file)$counts, keep.rownames= T)
           colnames(x)[2] <- "counts"
           x[]
         }, (.c)]
         DF <- dcast(.c, rn~prefix_output, value.var = "counts")
         setcolorder(DF, c("rn", rownames(sampleTable)))
         DF <- data.frame(DF, row.names = "rn")
         DF <- DF[rowSums(DF)>10,]
         dds <- DESeqDataSetFromMatrix(countData= DF, 
                                       colData= sampleTable, 
                                       design= ~rep+cdition)
         dds <- DESeq(dds)
         saveRDS(dds, dds_file)
         print(paste0(dds_file, "--> DONE!"))
       }else
         print(paste0(dds_file, "--> ALREADY EXISTS!"))
}, DESeq2_object]

#----------------------------------------------------------#
# Differential expression
#----------------------------------------------------------#
dir.create("db/FC_tables/", showWarnings = F)
dds <- meta[type=="RNA"]
dds <- dds[, .(DESeq2_object= unlist(tstrsplit(DESeq2_object, ";"))), setdiff(colnames(dds), "DESeq2_object")]
dds[, 
     {
       dds_file <- paste0("db/dds/", DESeq2_object, "_dds.rds")
       c_dds <- readRDS(dds_file)
       cmb <- as.character(c_dds$cdition)
       cmb <- unique(CJ(V1= cmb, V2= cmb))
       cmb <- cmb[V1!=V2]
       # cmb[, DESeq2_object:= DESeq2_object]
       cmb[, {
         FC_file <- paste0("db/FC_tables/", DESeq2_object, "_", V1, "_vs_", V2, ".txt")
         if(!file.exists(FC_file))
         {
           .c <- as.data.frame(lfcShrink(c_dds,
                                         type= "ashr",
                                         contrast= c("cdition", V1, V2)))
           fwrite(.c, 
                  FC_file, 
                  col.names = T, 
                  row.names = T, 
                  sep ="\t", 
                  quote= F)
           print(paste(FC_file, "--> DONE!"))
         }
         else 
           print(paste(FC_file, "--> ALREADY EXISTS!"))
       }, (cmb)]
}, DESeq2_object]


#----------------------------------------------------------#
# FPKMs
#----------------------------------------------------------#
dir.create("db/fpkms/", showWarnings = F)

dds <- meta[type=="RNA"]
dds <- dds[, .(DESeq2_object= unlist(tstrsplit(DESeq2_object, ";"))), setdiff(colnames(dds), "DESeq2_object")]

glength <- as.data.table(readRDS("db/counts/RNA_epiCancer/Ez18_1_counts.rds")$annotation)
glength <- glength[, .(GeneID, Length)]
setkeyv(glength, "GeneID")

dds[, 
     {
       fpkm_file <- paste0("db/fpkms/", DESeq2_object, ".txt")
       if(!file.exists(fpkm_file))
       {
         dds_file <- paste0("db/dds/", DESeq2_object, "_dds.rds")
         .c <- readRDS(dds_file)
         mcols(.c)$basepairs <- glength[names(.c@rowRanges), Length]
         res <- as.data.table(fpkm(.c), keep.rownames = "FBgn")
         fwrite(res, 
                fpkm_file, 
                col.names = T, 
                row.names = F, 
                sep= "\t", 
                quote = F, 
                na = NA)
         print(paste0(fpkm_file, "--> DONE!"))
       }else
         print(paste0(fpkm_file, "--> ALREADY EXISTS!"))
}, DESeq2_object]

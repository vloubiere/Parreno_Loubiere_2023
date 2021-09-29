setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

meta <- read_xlsx("Rdata/raw_metadata_final.xlsx")
meta <- as.data.table(meta)[type=="RNA"]

#----------------------------------------------------------#
# Build dm6 index
#----------------------------------------------------------#
if(!file.exists("/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index.log"))
{
  ref <- "/mnt/d/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
  buildindex(basename= "/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", reference= ref)
}

#----------------------------------------------------------#
# Alignment
#----------------------------------------------------------#
# Retrieve fq files starting with local, then external drive (in case where both exist, local will take over)
meta[, fq:= list.files(path = paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/fastq/", project),
                       pattern = fq_basename, 
                       recursive = T, 
                       full.names = T), .(project, fq_basename)]
meta[is.na(fq), fq:= list.files(path = paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/fastq/", project),
                                pattern = fq_basename, 
                                recursive = T, 
                                full.names = T), .(project, fq_basename)]
# Retrieve already created bam files
meta[, bam:= list.files(paste0("/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/", project), 
                        paste0(prefix_output, ".*.bam$"), 
                        recursive = T, 
                        full.names = T), .(project, prefix_output)]
# Generate names for the ones not created yet
meta[is.na(bam), 
     bam:= paste0("/mnt/d/_R_data/projects/epigenetic_cancer/db/bam/", project, "/", prefix_output, ".bam"), 
     .(project, prefix_output)]
# Create folders if not exist
meta[!file.exists(bam), 
     if(.N>0)
       dir.create(unique(dirname(bam)), showWarnings = F), project]
# Align
meta[!file.exists(bam), {
  if(.N>0)
  {
    print(paste0("START...", bam))
    subjunc(index= "/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index",
            readfile1= grep("_1.fq.gz$", fq, value = T),
            readfile2= if(layout=="PAIRED") grep("_2.fq.gz$", fq, value = T) else NULL,
            maxMismatches = 6,
            nthreads = 10,
            unique = T,
            output_file= bam)
    print("--> DONE!")
  }
}, .(bam, layout)]

#----------------------------------------------------------#
# Counts
#----------------------------------------------------------#
dir.create("db/counts/", showWarnings = F)
# Retrieve files for which counting has been already done
meta[, counts_file:= list.files(paste0("db/counts/", project),
                                paste0(prefix_output, "_counts.rds$"), 
                                full.names = T,
                                recursive = T), .(project, prefix_output)]
# Generate names for the ones not created yet
meta[is.na(counts_file), 
     counts_file:= paste0("db/counts/", project, "/", prefix_output, "_counts.rds"), 
     .(project, prefix_output)]
# Generate name for txt version
meta[!is.na(counts_file), counts_file_txt:= gsub(".rds$", "_only.txt", counts_file)]
# Create folders if not exist
meta[!file.exists(counts_file), 
     if((.N>0))
       dir.create(unique(dirname(counts_file)), showWarnings = F), project]
# Count
meta[, {
  if(.N>0)
  {
    # Count using rsubread
    if(!file.exists(counts_file))
    {
      .c <- featureCounts(bam,
                          annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf",
                          isGTFAnnotationFile = T,
                          isPairedEnd = layout=="PAIRED",
                          nthreads = 8)
      saveRDS(.c, counts_file)
    }else
      .c <- readRDS(counts_file)
    # Counts only text file
    if(!file.exists(counts_file_txt))
    {
      .c <- as.data.table(.c$counts, keep.rownames= T)
      names(.c) <- c("FBgn", "counts")
      fwrite(.c, 
             counts_file_txt, 
             sep= "\t", 
             quote= F)
      
    }
    print(paste(counts_file, "--> DONE!"))
  }
}, .(counts_file, counts_file_txt, bam, layout)]
# Text file version

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
dir.create("db/dds/", showWarnings = F)
meta[, dds_file:= paste0("db/dds/", DESeq2_object, ".rds"), DESeq2_object]
# Make list usable FBgns
genes <- import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
genes <- as.data.table(genes)
list_FBgn <- unique(genes[type %in% c("mRNA", "ncRNA"), gene_id])
meta[!file.exists(dds_file), 
     {
       if(.N>0)
       {
         sampleTable <- data.frame(unique(.SD[, .(cdition, rep, prefix_output, counts_file_txt)]), 
                                   row.names = "prefix_output")
         # Checks
         print(unique(project))
         print(unique(tissue))
         print(unique(method))
         print(unique(prep))
         print(sampleTable)
         N_reps <- any(sapply(split(sampleTable, sampleTable$cdition), nrow)<2)
         if(N_reps)
           print("Not enough replicates --> SKIPPED\n") else
           {
             DF <- lapply(sampleTable$counts_file_txt, fread)
             names(DF) <- rownames(sampleTable)
             DF <- rbindlist(DF, idcol = T)
             DF <- DF[FBgn %in% list_FBgn]
             DF <- dcast(DF, 
                         FBgn~.id, 
                         value.var = "counts")
             setcolorder(DF, c("FBgn", rownames(sampleTable)))
             DF <- data.frame(DF, row.names = "FBgn")
             DF <- DF[rowSums(DF)>ifelse(ncol(DF)<10, 10, ncol(DF)),]
             dds <- DESeqDataSetFromMatrix(countData= DF,
                                           colData= sampleTable,
                                           design= ~rep+cdition)
             dds <- DESeq(dds)
             saveRDS(dds, dds_file)
             print(paste0(dds_file, "--> DONE!"))
           }
       }
     }, dds_file]

#----------------------------------------------------------#
# Differential expression
#----------------------------------------------------------#
dir.create("db/FC_tables/", showWarnings = F)
meta[file.exists(dds_file),
     {
       c_dds <- readRDS(dds_file)
       cmb <- as.character(c_dds$cdition)
       cmb <- unique(CJ(V1= cmb, V2= cmb))
       cmb <- cmb[V1!=V2]
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
     }, .(DESeq2_object, dds_file)]

#----------------------------------------------------------#
# FPKMs
#----------------------------------------------------------#
dir.create("db/fpkms/", showWarnings = F)
glength <- as.data.table(readRDS("db/counts/RNA_epiCancer/Ez18_1_counts.rds")$annotation)
glength <- glength[, .(GeneID, Length)]
setkeyv(glength, "GeneID")
meta[, fpkm_file:= paste0("db/fpkms/", DESeq2_object, ".txt"), DESeq2_object]
meta[file.exists(dds_file) & !file.exists(fpkm_file), 
     {
       if(.N>0)
       {
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
       }
     }, .(fpkm_file, dds_file, DESeq2_object)]

#----------------------------------------------------------#
# BW files
#----------------------------------------------------------#
dir.create("db/bw", 
           showWarnings = F)
meta[project=="RNA_epiCancer" & prep=="Parreno" & rep==1 & method=="handDissect" & tissue=="ED", 
     bw_file:= paste0("db/bw/", DESeq2_object, "_", cdition, ".bw"), .(DESeq2_object, cdition)]

meta[!is.na(bw_file) & !file.exists(bw_file),
     {
       if(.N>0)
       {
         print("START...")
         .c <- fread(cmd= paste0("/usr/local/bin/samtools view -@ 10 ", bam), 
                     sep = "\t", 
                     fill = T, 
                     header = F, 
                     sel = c(2,3,4,5,9))
         .c <- .c[V2 %in% c(99, 163) & V5>=30, .(seqnames= V3,
                                                 start= V4,
                                                 end= V4+V9)]
         gr <- GRanges(.c)
         cov <- GenomicRanges::coverage(gr)/length(gr)*1e6
         export.bw(GRanges(cov), 
                   con= bw_file)
         print(paste(bw_file, "--> DONE!"))
       }
     }, .(bam, bw_file)]

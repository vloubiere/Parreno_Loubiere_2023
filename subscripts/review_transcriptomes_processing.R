setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(readxl)
require(Rsubread)

# Import metadata ----
meta <- as.data.table(read_xlsx("Rdata/metadata_RNA_review.xlsx"))
meta <- meta[!is.na(DESeq2_object)]

# Create output files ----
meta[, bam:= paste0("/scratch/stark/vloubiere/epicancer/",
                    DESeq2_object, "_", cdition, "_rep", rep, ".bam"), .(DESeq2_object, cdition, rep)]
meta[, read_counts:= paste0("db/counts/RNA/", DESeq2_object, "_read_counts.rds"), DESeq2_object]
meta[, dds_file:= paste0("db/dds/RNA/", DESeq2_object, "_dds.rds"), DESeq2_object]

# Retrieve fq files ----
meta[, c("fq1", "fq2"):= lapply(.SD, function(x){
  list.files("db/fq/RNA/", 
             paste0("^", x, "$"), 
             recursive= T,
             full.names = T)
}), by= .(fq1, fq2, project), .SDcols= c("fq1", "fq2")]

# BUILD dm6 index  ----
if(!file.exists("/groups/stark/vloubiere/genomes/Drosophila_melanogaster/subreadr_dm6/subreadr_dm6_index.log"))
{
  buildindex(basename= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/subreadr_dm6/subreadr_dm6_index",
             reference= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/UCSC/dm6/Sequence/WholeGenomeFasta/genome.fa")
}

# ALIGNMENT ----
meta[, {
  if(!file.exists(bam))
  {
    print(paste0("START...", bam))
    if(.N>1) # Merge fastq files when re-sequencing
    {
      tmp <- tempfile(tmpdir = "/scratch/stark/vloubiere/lauriane/", fileext = "_1.gz")
      cmd <- paste(c("cat", paste0("'", normalizePath(fq1), "'"), ">", tmp), collapse = " ")
      system(cmd)
      fq1 <- tmp
      tmp <- tempfile(tmpdir = "/scratch/stark/vloubiere/lauriane/", fileext = "_2.gz")
      cmd <- paste(c("cat", paste0("'", normalizePath(fq2), "'"), ">", tmp), collapse = " ")
      system(cmd)
      fq2 <- tmp
    }
    align(index= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/subreadr_dm6/subreadr_dm6_index",
          readfile1= fq1,
          readfile2= fq2,
          maxMismatches = 6,
          nthreads = 10,
          unique = T,
          output_file= bam)
  }
  print(paste(bam, "--> DONE!"))
}, bam]

# Counts ----
meta[, {
  if(!file.exists(read_counts))
  {
    .c <- featureCounts(unique(bam), # count reads
                        annot.ext= "/groups/stark/vloubiere/genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf",
                        isGTFAnnotationFile = T,
                        isPairedEnd = T,
                        nthreads = 8)
    saveRDS(.c, read_counts)
  }
}, read_counts]

# DESeq2 ----
FC <- meta[, {
  if(!file.exists(dds_file))
  {
    # SampleTable and DF
    .c <- readRDS(read_counts)
    DF <- as.data.frame(.c[[1]])
    DF <- DF[order(rownames(DF)),]
    colnames(DF) <- unique(paste0(cdition, "__rep", rep))
    sampleTable <- as.data.frame(setNames(tstrsplit(names(DF), "__"),
                                          c("cdition", "rep")),
                                 row.names = names(DF))
    # Filter low read genes
    DF <- DF[rowSums(DF)>10, ] # filter low read genes + N
    # dds
    dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                          colData= sampleTable,
                                          design= ~rep+cdition)
    dds <- DESeq2::DESeq(dds)
    # Add gene length to compute fpkms
    glength <- as.data.table(.c$annotation, key= "GeneID")[, .(GeneID, Length)]
    GenomicRanges::mcols(dds)$basepairs <- glength[names(dds@rowRanges), Length]
    # SAVE
    saveRDS(dds, dds_file)
  }else
    dds <- readRDS(dds_file)
  CJ(cdition= as.character(dds$cdition),
     control= as.character(dds$cdition), 
     unique = T)
}, .(read_counts, dds_file, DESeq2_object)]

# FC tables ----
cditions_table <- read_xlsx("Rdata/RNA_dds_conditions_review.xlsx")
cditions_table <- as.data.table(cditions_table)
FC <- FC[cditions_table, on= c("cdition", "control"), nomatch= NULL]
FC[, FC_file:= paste0("db/FC_tables/RNA/", DESeq2_object, "_", cdition, "_vs_", control, ".txt")]
FC[, {
  dds <- readRDS(dds_file)
  .SD[, {
    if(!file.exists(FC_file))
    {
      # Compute FC tables
      .c <- DESeq2::lfcShrink(dds,
                              type= "ashr",
                              contrast= c("cdition", cdition, control))
      .c <- as.data.table(as.data.frame(.c), 
                          keep.rownames = "FBgn")
      .c[, diff:= fcase(padj<0.05 & log2FoldChange>1, "up",
                        padj<0.05 & log2FoldChange<(-1), "down",
                        default= "unaffected")]
      fwrite(.c,
             FC_file, 
             col.names = T,
             sep= "\t")
    }
    print("DONE")
  }, .(cdition, control, FC_file)]
  print("DONE")
}, dds_file]

# Add to meta ----
meta <- merge(meta, 
              FC[, .(FC_file= paste0(FC_file, collapse = ",")), cdition],
              by= "cdition", 
              all.x= T)

# SAVE ----
fwrite(meta, 
       "Rdata/processed_metadata_RNA_review.txt",
       na = NA)

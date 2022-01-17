setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

meta <- as.data.table(read_xlsx("Rdata/metadata_RNA.xlsx"))
meta <- meta[DESeq2_object %in% c("epiCancer_ED_allograft_RNA_gDNA", 
                                  "epiCancer_ED_RNA_CUTNRUN",
                                  "RNA_development_Public",
                                  "RNA_mutants_SA2020_ED_Martinez",
                                  "RNA_Paro_2018_Paro",
                                  "RNA_phRNAi_SA2020_ED_Loubiere")]
# Retrieve fq files starting with local, then external drive (in case where both exist, local will take over)
hdir <- "/mnt/f/_R_data/projects/epigenetic_cancer/db/"
meta[, c("fq1", "fq2"):= lapply(.SD, function(x){
  list.files(paste0(hdir, "fastq/", project), 
             paste0("^", x, "$"), 
             recursive= T,
             full.names = T)
}), by= .(fq1, fq2, project), .SDcols= c("fq1", "fq2")]

#----------------------------------------------------------#
# Output files
#----------------------------------------------------------#
meta[, bam:= paste0(hdir, "bam/", project, "/", prefix_output, ".bam"), .(project, prefix_output)]
meta[, dds_file:= paste0("db/dds/", DESeq2_object, "_dds.rds"), DESeq2_object]
fwrite(meta, "Rdata/processed_metadata_RNA.txt")

#----------------------------------------------------------#
# ALIGNMENT
#----------------------------------------------------------#
# Build dm6 index
if(!file.exists("/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index.log"))
{
  ref <- "/mnt/d/_R_data/genomes/dm6/Sequence/WholeGenomeFasta/genome.fa"
  buildindex(basename= "/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index", reference= ref)
}

# Align
meta[, {
  if(!file.exists(bam))
  {
    print(paste0("START...", bam))
    subjunc(index= "/mnt/d/_R_data/genomes/dm6/subreadr_index/subreadr_dm6_index",
            readfile1= fq1,
            readfile2= ifelse(is.na(fq2), NULL, fq2),
            maxMismatches = 6,
            nthreads = 10,
            unique = T,
            output_file= bam)
    
  }
  print(paste(bam, "--> DONE!"))
}, .(bam, fq1, fq2)]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
# DESEQ
meta[, {
  if(!file.exists(dds_file))
  {
    # SampleTable and DF
    sampleTable <- .SD[, .(check= .N>1, rep= paste0("rep", rep), bam, paired= !is.na(fq2)), cdition][(check)]
    .c <- featureCounts(sampleTable$bam, # count reads
                        annot.ext= "../../genomes/dm6/dmel-all-r6.36.gtf",
                        isGTFAnnotationFile = T,
                        isPairedEnd = sampleTable$paired,
                        nthreads = 10)
    DF <- as.data.frame(.c[[1]])
    sampleTable <- data.frame(sampleTable[, .(cdition, rep)], 
                              row.names = basename(sampleTable$bam))
    # Filter low read genes
    DF <- DF[rowSums(DF)>10, ] # filter low read genes + N
    # dds
    dds <- DESeq2::DESeqDataSetFromMatrix(countData= DF,
                                          colData= sampleTable,
                                          design= ~rep+cdition)
    dds <- DESeq2::DESeq(dds)
    # Add gene length to compute fpkms
    glength <- as.data.table(.c$annotation, key= "GeneID")[, .(GeneID, Length)]
    mcols(dds)$basepairs <- glength[names(dds@rowRanges), Length]
    # SAVE
    saveRDS(dds, dds_file)
  }else{
    # FC tables
    dds <- readRDS(dds_file)
    cmb <- CJ(dds$cdition, dds$cdition, unique = T)
    cmb <- cmb[V1!=V2]
    cmb[, FC_file:= paste0("db/FC_tables/", DESeq2_object, "_", V1, "_vs_", V2, ".txt")]
    cmb[, {
      if(!file.exists(FC_file))
      {
        .c <- DESeq2::lfcShrink(dds,
                                type= "ashr",
                                contrast= c("cdition",
                                            as.character(V1),
                                            as.character(V2)))
        .c <- suppressWarnings(as.data.table(.c, keep.rownames= "FBgn"))
        fwrite(.c,
               FC_file,
               col.names = T,
               row.names = T,
               sep ="\t",
               quote= F)
      }
    }, .(V1, V2, FC_file)]
  }
  print("DONE")
}, .(DESeq2_object, dds_file)]

#----------------------------------------------------------#
# BW files
#----------------------------------------------------------#
# dir.create("db/bw", 
#            showWarnings = F)
# meta[project=="RNA_epiCancer" & prep=="Parreno" & rep==1 & method=="handDissect" & tissue=="ED", 
#      bw_file:= paste0("db/bw/", DESeq2_object, "_", cdition, ".bw"), .(DESeq2_object, cdition)]
# 
# meta[!is.na(bw_file) & !file.exists(bw_file),
#      {
#        if(.N>0)
#        {
#          print("START...")
#          .c <- fread(cmd= paste0("/usr/local/bin/samtools view -@ 10 ", bam), 
#                      sep = "\t", 
#                      fill = T, 
#                      header = F, 
#                      sel = c(2,3,4,5,9))
#          .c <- .c[V2 %in% c(99, 163) & V5>=30, .(seqnames= V3,
#                                                  start= V4,
#                                                  end= V4+V9)]
#          gr <- GRanges(.c)
#          cov <- GenomicRanges::coverage(gr)/length(gr)*1e6
#          export.bw(GRanges(cov), 
#                    con= bw_file)
#          print(paste(bw_file, "--> DONE!"))
#        }
#      }, .(bam, bw_file)]

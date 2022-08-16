setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(GenomicRanges)
require(readxl)

meta <- as.data.table(read_xlsx("Rdata/metadata_RNA.xlsx"))
meta <- meta[DESeq2_object %in% c("epiCancer_ED_GFP+_system_RNA", 
                                  "epiCancer_ED_GFP-_system_RNA",
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
meta[, dds_file:= paste0("db/dds/RNA/", DESeq2_object, "_dds.rds"), DESeq2_object]

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
# BW files
#----------------------------------------------------------#
meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA", bw_file:= paste0("db/bw/RNA-Seq/", gsub(".bam", ".bw", basename(bam)))]
meta[!is.na(bw_file), {
  if(!file.exists(bw_file))
  {
    bed <- fread(cmd= paste0("bamToBed -i ", bam, " -bedpe"), 
                 fill= T, sep= "\t")
    bed <- bed[V1==V4 & V6>V2 & V8>=30, .(seqnames= V1, start= V2, end= V6)]
    bed <- GRanges(bed)
    cov <- coverage(bed)/length(bed)*1e6
    rtracklayer::export.bw(GRanges(cov), 
                           con= bw_file)
    print("done")
  }
}, bw_file]
meta[DESeq2_object=="epiCancer_ED_GFP-_system_RNA", bw_merge:= paste0("db/bw/RNA-Seq/", cdition, "_log2_merge.bw")]
meta[!is.na(bw_merge), {
  if(!file.exists(bw_merge))
    vl_bw_merge(bw_file, 
                output = bw_merge, 
                genome = "dm6", 
                scoreFUN = function(x) {
                  x <- x/length(unique(bw_file))
                  ifelse(x==0, x, log2(x))
                })
  print("done")
}, bw_merge]

#----------------------------------------------------------#
# DESeq2
#----------------------------------------------------------#
FC <- meta[, {
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
  }else
    dds <- readRDS(dds_file)
  FC <- CJ(as.character(dds$cdition), 
           as.character(dds$cdition), 
           unique = T)
  FC <- FC[V1!=V2]
}, .(DESeq2_object, dds_file)]

#----------------------------------------------------------#
# FC tables
#----------------------------------------------------------#
cditions_table <- read_xlsx("Rdata/RNA_dds_conditions.xlsx")
setkeyv(FC, c("V1", "V2"))
FC <- FC[as.data.table(cditions_table)]
FC[, FC_file:= paste0("db/FC_tables/RNA/", DESeq2_object, "_", V1, "_vs_", V2, ".txt"), .(DESeq2_object, V1, V2)]
FC[, {
  dds <- readRDS(dds_file)
  .SD[, {
    if(!file.exists(FC_file))
    {
      # Compute FC tables
      .c <- DESeq2::lfcShrink(dds,
                              type= "ashr",
                              contrast= c("cdition", V1, V2))
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
  }, .(V1, V2, FC_file)]
  print("DONE")
}, .(dds_file)]

meta[FC[V2!="RNA_PH29"], FC_file:= i.FC_file, on= c("DESeq2_object", "cdition==V1")]
meta[FC[V2=="RNA_PH29"], FC_file_PH29:= i.FC_file, on= c("DESeq2_object", "cdition==V1")]
fwrite(meta, 
       "Rdata/processed_metadata_RNA.txt",
       na = NA)


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

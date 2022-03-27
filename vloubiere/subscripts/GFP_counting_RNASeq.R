setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

meta <- fread("Rdata/processed_metadata_RNA.txt") 
meta <- meta[DESeq2_object %in% c("epiCancer_ED_RNA_CUTNRUN",
                                  "epiCancer_ED_allograft_RNA_gDNA")]

#----------------------------------------------------------#
# Output files
#----------------------------------------------------------#
meta[, bam:= paste0(switch(DESeq2_object,
                           "epiCancer_ED_RNA_CUTNRUN"= "db/bam/cutnrun_GFP/",
                           "epiCancer_ED_allograft_RNA_gDNA"= "db/bam/allograft_GFP/"), 
                    prefix_output, "_GFP.bam"), .(DESeq2_object, prefix_output)]

#----------------------------------------------------------#
# INDEXING
#----------------------------------------------------------#
# Build GFP idx (https://www.addgene.org/74749/sequences/)
if(!file.exists("/mnt/d/_R_data/genomes/nowGFP/GFP_index.log"))
{
  seq <- "gtgagcaagggcgagaagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgtccggcgagggcgagggcgatgccacctacggcaagatgtccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcaaaaccaccctgacctggggcatgcagtgcttcgcccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacgtccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcgtcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacgccatcagcggcaacgccaatatcaccgccgacaagcagaagaacggcatcaaggcctacttcacgatccgccacgacgtcgaggacggcagcgtgctgctcgccgaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcacccagtccaagcagagcaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatccctctcggcgcggacgagctgtacaag"
  seqinr::write.fasta(sequences = seq,
                      names= "GFP", 
                      as.string = T, 
                      file.out = "/mnt/d/_R_data/genomes/nowGFP/GFP.fa")
  ref <- "/mnt/d/_R_data/genomes/nowGFP/GFP.fa"
  buildindex(basename= "/mnt/d/_R_data/genomes/nowGFP/GFP_index", 
             reference= ref)
}
if(!file.exists("/mnt/d/_R_data/genomes/GFP_RFP/GFP_index.log"))
{
  sequences <- c(EGFP= "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTAA",
                 mRFP1= "ATGGCCTCCTCCGAGGACGTCATCAAGGAGTTCATGCGCTTCAAGGTGCGCATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGCGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCCAGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCACCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGATGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCCGAGGTCAAGACCACCTACATGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAAGACCGACATCAAGCTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAGCGCGCCGAGGGCCGCCACTCCACCGGCGCCTAA")
  seqinr::write.fasta(sequences = as.list(sequences),
                      names= names(sequences), 
                      as.string = T, 
                      file.out = "/mnt/d/_R_data/genomes/GFP_RFP/GFP.fa")
  ref <- "/mnt/d/_R_data/genomes/GFP_RFP/GFP.fa"
  buildindex(basename= "/mnt/d/_R_data/genomes/GFP_RFP/GFP_index", 
             reference= ref)
}

#----------------------------------------------------------#
# ALIGNMENT
#----------------------------------------------------------#
# Align
meta[, index:= switch(DESeq2_object,
                      "epiCancer_ED_RNA_CUTNRUN"= "/mnt/d/_R_data/genomes/nowGFP/GFP_index",
                      "epiCancer_ED_allograft_RNA_gDNA"= "/mnt/d/_R_data/genomes/GFP_RFP/GFP_index"), DESeq2_object]
meta[, {
  if(!file.exists(bam))
  {
    print(paste0("START...", bam))
    subjunc(index= index,
            readfile1= fq1,
            readfile2= ifelse(is.na(fq2), NULL, fq2),
            maxMismatches = 6,
            nthreads = 10,
            unique = T,
            output_file= bam)
  }
  print(paste(bam, "--> DONE!"))
}, .(bam, fq1, fq2, index)]

#----------------------------------------------------------#
# IMPORT counts GFP reads
#----------------------------------------------------------#
if(!file.exists("db/counts/GFP_RFP/GFP_RFP_counts_epiCancer_RNA.txt"))
{
  counts <- meta[, {
    fread(cmd= paste0("/usr/local/bin/samtools view -@ 9 -b -q 30 ",
                      bam, " | bedtools bamtobed -i stdin"))
  }, .(DESeq2_object, dds_file, cdition, bam)]
  fwrite(counts, 
         "db/counts/GFP_RFP/GFP_RFP_counts_epiCancer_RNA.txt", 
         sep= "\t",
         na = NA)
}
counts <- fread("db/counts/GFP_RFP/GFP_RFP_counts_epiCancer_RNA.txt")
counts <- counts[, .N, .(cdition, dds_file, DESeq2_object, bam, V1)]

#----------------------------------------------------------#
# Add dds sizeFactors
#----------------------------------------------------------#
sizeFactors <- counts[, as.data.table(readRDS(dds_file)$sizeFactor, keep.rownames = T), dds_file]
sizeFactors[, V1:= gsub(".bam", "_GFP.bam", V1)]
counts[, bam_basename:= basename(bam)]
counts[sizeFactors, sizeFactor:= i.V2, on= c("bam_basename==V1", "dds_file")]
counts <- na.omit(counts[, .(norm_counts= .(N/sizeFactor),
                             mean= mean(N/sizeFactor)), .(DESeq2_object, cdition, V1)])

pdf("pdf/RNA/GFP_norm_counts.pdf", 
    width = 3,
    height= 9)
par(mfrow= c(3,1),
    mar= c(10,6,4,3))
# CUTNRUN genotype
bar <- barplot(counts[DESeq2_object=="epiCancer_ED_RNA_CUTNRUN", mean],
               las= 2,
               names.arg = counts[DESeq2_object=="epiCancer_ED_RNA_CUTNRUN", cdition],
               ylim= c(0, 27500),
               col= "limegreen")
mtext("GFP CUTNRUN genotype", line= 2)
mtext("GFP normalized counts", 2, line = 4)
points(rep(bar, each= 3),
       unlist(counts[DESeq2_object=="epiCancer_ED_RNA_CUTNRUN", norm_counts]))
# CUTNRUN genotype RFP
bar <- barplot(counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="mRFP1", mean],
               las= 2,
               names.arg = counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="mRFP1", cdition],
               ylim= c(0, 100000),
               col= "tomato")
mtext("mRFP1 ALLOGRAFT genotype", line= 2)
mtext("mRFP1 normalized counts", 2, line = 4)
reps <- counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="mRFP1", norm_counts]
points(rep(bar, lengths(reps)),
       unlist(counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="mRFP1", norm_counts]))
# CUTNRUN genotype GFP
bar <- barplot(counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="EGFP", mean],
               las= 2,
               names.arg = counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="EGFP", cdition],
               ylim= c(0, 10000),
               col= "limegreen")
mtext("EGFP ALLOGRAFT genotype", line= 2)
mtext("EGFP normalized counts", 2, line = 4)
reps <- counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="EGFP", norm_counts]
points(rep(bar, lengths(reps)),
       unlist(counts[DESeq2_object=="epiCancer_ED_allograft_RNA_gDNA" & V1=="EGFP", norm_counts]))

dev.off()
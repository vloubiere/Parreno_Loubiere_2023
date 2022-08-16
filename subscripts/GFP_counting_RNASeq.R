setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)

meta <- fread("Rdata/processed_metadata_RNA.txt") 
meta <- meta[DESeq2_object %in% c("epiCancer_ED_GFP-_system_RNA",
                                  "epiCancer_ED_GFP+_system_RNA")]

#----------------------------------------------------------#
# Output files
#----------------------------------------------------------#
meta[, bam:= paste0(switch(DESeq2_object,
                           "epiCancer_ED_GFP-_system_RNA"= "/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/cutnrun_GFP/",
                           "epiCancer_ED_GFP+_system_RNA"= "/mnt/f/_R_data/projects/epigenetic_cancer/db/bam/allograft_GFP/"), 
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
                      "epiCancer_ED_GFP-_system_RNA"= "/mnt/d/_R_data/genomes/nowGFP/GFP_index",
                      "epiCancer_ED_GFP+_system_RNA"= "/mnt/d/_R_data/genomes/GFP_RFP/GFP_index"), DESeq2_object]
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
    fread(cmd= paste0("/usr/bin/samtools view -@ 9 -b -q 30 ",
                      bam, " | bedtools bamtobed -i stdin"))
  }, .(DESeq2_object, dds_file, cdition, bam)]
  fwrite(counts, 
         "db/counts/GFP_RFP/GFP_RFP_counts_epiCancer_RNA.txt", 
         sep= "\t",
         na = NA)
}
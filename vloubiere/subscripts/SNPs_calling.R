setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)
require(Rsubread)

#----------------------------------#
# Merge fastq file
#----------------------------------#
dat <- fread("Rdata/raw_metadata.txt")[project=="DNA_epiCancer"]
dat[grepl("_1.fq.gz$", fq_file), merged_file:= ifelse(.N>1, paste0("db/fastq/", dirname(fq_file), "/", cdition, "_merged_1.fq.gz"), fq_file), cdition]
dat[grepl("_2.fq.gz$", fq_file), merged_file:= ifelse(.N>1, paste0("db/fastq/", dirname(fq_file), "/", cdition, "_merged_2.fq.gz"), fq_file), cdition]
dat[, {
  if(!file.exists(merged_file) & .N>1)
  {
    cmd <- paste("cat", paste(fq, collapse= " "), ">", merged_file)
    system(vl_bash_wrap_windows(cmd))
    print(cmd)
  }
}, .(cdition, merged_file)]

#----------------------------------#
# Alignment
#----------------------------------#
dat <- fread("Rdata/raw_metadata.txt")[project %in% c("DNA_epiCancer", "DNA_tano")]
dat[, sam:= paste0("db/sam/", project, "/", cdition, ".sam")]
dat[, {
  if(!file.exists(sam))
  {
    align(index= "../../genomes/dm6/subreadr_index/subreadr_dm6_index",
          readfile1= paste0("db/fastq/", grep("_1.fq.gz$", fq_file, value = T)),
          readfile2 = paste0("db/fastq/", grep("_2.fq.gz$", fq_file, value = T)),
          type = "dna",
          output_format = "SAM",
          output_file = sam,
          maxMismatches = 6,
          unique = T,
          indels = 5,
          complexIndels = T,
          nTrim5 = 3,
          nthreads = 10,
          # detect structural variants
          detectSV = T)
  }
}, .(cdition, sam)]

#----------------------------------#
# Remove duplicates
#----------------------------------#
dat[, coll_sam:= gsub(".sam$", "_uniq.sam", sam)]
dat <- dat[project=="DNA_tano" & cdition %in% as.factor(x = c("LGL_Donor", "LGL_T10"))]
dat[, {
  if(!file.exists(coll_sam))
  {
    removeDupReads(sam,
                   threshold = 1,
                   coll_sam,
                   outputFormat = "SAM")
  }
}, .(sam, coll_sam)]


#----------------------------------#
# SNP ref genome
#----------------------------------#
dat[, vcf := paste0("db/vcf/", gsub(".sam", "_vs_refgenome.vcf", basename(sam))), sam]
dat[, {
  if(!file.exists(vcf))
  {
    Rsubread::exactSNP(readFile = sam, 
                       refGenomeFile = "../../genomes/dm6/Sequence/WholeGenomeFasta/genome.fa", 
                       outputFile = vcf,
                       nthreads = 10)
  }
  print("")
}, .(sam, vcf)] 

# combinations <- data.table(cdition= c("ph18", "ph18", "ph29", "ph29", "ph29_t5", "phd11", "phd11_t8"), 
#                            annotation= c("refgenome", "ph18", "refgenome", "ph18", "refgenome", "refgenome", "refgenome"))
# combinations[, sam:= paste0("db/sam/DNA_epicancer_VL/", cdition, ".sam")]
# combinations[, vcf:= paste0("db/vcf/", cdition, "_vs_", annotation, ".vcf")]
# combinations[annotation!="refgenome", annotation_vcf:= paste0("db/vcf/", annotation, "_vs_refgenome.vcf")]
# combinations[, {
#   if(!file.exists(vcf))
#   {
#     if(is.na(annotation_vcf))
#     {
#       Rsubread::exactSNP(readFile = sam, 
#                          refGenomeFile = "../../genomes/dm6/Sequence/WholeGenomeFasta/genome.fa", 
#                          outputFile = vcf,
#                          nthreads = 10)
#     }else{
#       Rsubread::exactSNP(readFile = sam, 
#                          refGenomeFile = "../../genomes/dm6/Sequence/WholeGenomeFasta/genome.fa", 
#                          SNPAnnotationFile = annotation_vcf,
#                          outputFile = vcf,
#                          nthreads = 10)
#     }
#   }
#   print("")
# }, (combinations)] 

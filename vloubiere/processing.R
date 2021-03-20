setwd("D:/_R_data/projects/epigenetic_cancer/")
sapply(list.files("/_R_data/functions/", ".R$", full.names = T), source)
require(Rsubread)
require(DESeq2)
require(data.table)
require(rtracklayer)
require(kohonen)
require(pheatmap)
require(ontologyIndex)
require(digest)

if(1==0)
{
  source("git_epiCancer/vloubiere/subscripts/alignment_statistics.R")
  source("git_epiCancer/vloubiere/subscripts/clustering.R")
  source("git_epiCancer/vloubiere/subscripts/compute_transcripts_counts.R")
  source("git_epiCancer/vloubiere/subscripts/DESEQ2_analysis.R")
  source("git_epiCancer/vloubiere/subscripts/dm6_index_Rsubread.R")
  source("git_epiCancer/vloubiere/subscripts/Gene_ontologies.R")
  source("git_epiCancer/vloubiere/subscripts/PcG_binding_binary.R")
  source("git_epiCancer/vloubiere/subscripts/Processing_ATAC_seq.R")
  source("git_epiCancer/vloubiere/subscripts/Public_ChIP_quantif_clusters.R")
  source("git_epiCancer/vloubiere/subscripts/Public_transcriptomes_quantif_clusters.R")
  source("git_epiCancer/vloubiere/subscripts/String_networks.R")
  source("git_epiCancer/vloubiere/subscripts/Targeted_heatmaps.R")
  source("git_epiCancer/vloubiere/subscripts/transcriptomes_alignment.R")
  source("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R")
}

# Splicing
if(1==0)
{
  source("git_epiCancer/vloubiere/subscripts/polarity_GO_JP.R")
}

  
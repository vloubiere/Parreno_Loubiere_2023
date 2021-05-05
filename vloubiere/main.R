setwd("D:/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(Rsubread)
require(DESeq2)
require(data.table)
require(rtracklayer)
require(kohonen)
require(pheatmap)
require(ontologyIndex)
require(digest)
require(gridExtra)
require(readxl)
require(STRINGdb)
require(TFBSTools)
require(BSgenome.Dmelanogaster.UCSC.dm6)
require(pdftools)
require(diagram)

#### Transcritpomes raw data processing
if(F)
{
  file.edit("git_epiCancer/vloubiere/subscripts/dm6_index_Rsubread.R")
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/download_modEncode_RNASeq.R") # Not working
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/download_Paro_data.R")
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/transcriptomes_processing.R")
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/MA_plots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/alignment_statistics.R")
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R")
}

##### Compute available data
if(F)
{
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_ATAC_seq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_K9me3_seq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/REs_gene_assignment.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Compute_available_data.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Compute_motifs.R")
  file.edit("git_epiCancer/vloubiere/subscripts/tissue_specificity_score.R")
}

#### Clustering transcriptomes
if(F)
{
  file.edit("git_epiCancer/vloubiere/subscripts/clustering.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Targeted_heatmaps.R")  
  file.edit("git_epiCancer/vloubiere/subscripts/Gene_ontologies.R")
  file.edit("git_epiCancer/vloubiere/subscripts/String_networks.R")
  
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_heatmaps.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_transriptomes.R")
}

#### Analysis transplantations
if(F)
{
  file.edit("git_epiCancer/vloubiere/subscripts/transplantation_experiments.R")
}

#### Further analyses
if(F)
{
  file.edit("git_epiCancer/vloubiere/subscripts/motif_wt_activity.R") # Aim was to classify as rep/act. Not very conclusive
  file.edit("git_epiCancer/vloubiere/subscripts/compare_Paro_data.R")
}

# Splicing for JP
if(1==0)
{
  file.edit("git_epiCancer/vloubiere/subscripts/polarity_GO_JP.R")
  file.edit("git_epiCancer/vloubiere/subscripts/motifs_enrichment_AMaria.R")
}

  
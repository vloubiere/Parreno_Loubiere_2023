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

#### Transcritpomes raw data processing
if(F)
{
  source("git_epiCancer/vloubiere/subscripts/dm6_index_Rsubread.R")
  source("git_epiCancer/vloubiere/subscripts/transcriptomes_alignment.R")
  source("git_epiCancer/vloubiere/subscripts/alignment_statistics.R")
  source("git_epiCancer/vloubiere/subscripts/compute_transcripts_counts.R")
  source("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R")
  source("git_epiCancer/vloubiere/subscripts/DESEQ2_analysis.R")
}

##### Compute available data
if(F)
{
  source("git_epiCancer/vloubiere/subscripts/Processing_ATAC_seq.R")
  source("git_epiCancer/vloubiere/subscripts/REs_gene_assignment.R")
  source("git_epiCancer/vloubiere/subscripts/Compute_available_data.R")
}

#### Clustering transcriptomes
if(F)
{
  source("git_epiCancer/vloubiere/subscripts/clustering.R")
  source("git_epiCancer/vloubiere/subscripts/Targeted_heatmaps.R")  
  source("git_epiCancer/vloubiere/subscripts/Gene_ontologies.R")
  source("git_epiCancer/vloubiere/subscripts/String_networks.R")
  
  source("git_epiCancer/vloubiere/subscripts/clustering_compute_additional_data.R")
  source("git_epiCancer/vloubiere/subscripts/clustering_heatmaps.R")
}


#### Clusters analysis
# detailed heatmaps

# Gene Ontologies
# STRING networks

# Add other data
source("git_epiCancer/vloubiere/subscripts/Processing_ATAC_seq.R")
source("git_epiCancer/vloubiere/subscripts/PcG_binding_gene_assignment.R")
source("git_epiCancer/vloubiere/subscripts/Public_ChIP_quantif_clusters.R")
source("git_epiCancer/vloubiere/subscripts/Public_transcriptomes_quantif_clusters.R")
# Plot overlay external data and clusters
file.edit("git_epiCancer/vloubiere/subscripts/plot_overlay_available_data_clusters.R")



# Splicing for JP
if(1==0)
{
  source("git_epiCancer/vloubiere/subscripts/polarity_GO_JP.R")
}

  
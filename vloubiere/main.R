setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

if(F)
{
  #-------------------------------------------------------------------------#
  # Available data
  #-------------------------------------------------------------------------#
  #### Process public data and identify/classify ED REs ####
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_available_data.R")
  file.edit("git_epiCancer/vloubiere/subscripts/ATAC_Seq_peak_calling.R")
  file.edit("git_epiCancer/vloubiere/subscripts/REs_assignment_and_features.R")
  
  #-------------------------------------------------------------------------#
  # RNA analyses
  #-------------------------------------------------------------------------#
  #### Process RNA-Seq data ####
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/alignment_statistics.R")
  file.edit("git_epiCancer/vloubiere/subscripts/generate_final_FC_table.R")
  
  #### Bulk analyses ####
  file.edit("git_epiCancer/vloubiere/subscripts/MA_plots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R")
  file.edit("git_epiCancer/vloubiere/subscripts/alluvial_plot_affected_genes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_transriptomes.R")
  
  #### Clustering ####
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_transcriptomes.R")
  
  #-------------------------------------------------------------------------#
  # Cut N run
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_QC.R")
  file.edit("git_epiCancer/vloubiere/subscripts/compare_vl_gonza_cutnrun_tracks.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_peak_calling.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_differential_analysis.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_diff_regions_analysis.R")
  file.edit("git_epiCancer/vloubiere/subscripts/GO_K27_cluster_genes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Screenshots_K27_clusters.R")
  file.edit("git_epiCancer/vloubiere/subscripts/K27_clusters_intersections.R")
  
  
  #-------------------------------------------------------------------------#
  # DNA analyses
  #-------------------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------------------#
  # MISCELLANEOUS / Casual
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/dac_locus_Giacomo.R")
} 



  
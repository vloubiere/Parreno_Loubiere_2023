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
  file.edit("git_epiCancer/vloubiere/subscripts/generate_final_FC_tables.R")
  
  #### Bulk analyses ####
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R")
  file.edit("git_epiCancer/vloubiere/subscripts/MA_plots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/alluvial_plot_affected_genes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_transriptomes.R")
  
  #### Clustering ####
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_allograft_transcriptomes_PH29_PHD11.R")
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_cutnrun_transcriptomes_PH29_PHD11.R")
  
  #-------------------------------------------------------------------------#
  # Cut N run
  #-------------------------------------------------------------------------#
  # Processing and QC
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_QC.R")
  file.edit("git_epiCancer/vloubiere/subscripts/compare_vl_gonza_cutnrun_tracks.R")
  
  # Analyses
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_peak_calling.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_K27_changes_volcano_plot.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_clustering_PH29_PHD11.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_clustering_PHD9_PHD11.R")
  
  #-------------------------------------------------------------------------#
  # Overlap cutnrun / transriptomes
  #-------------------------------------------------------------------------#
  
  
  #-------------------------------------------------------------------------#
  # DNA analyses
  #-------------------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------------------#
  # MISCELLANEOUS / Casual
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/dac_locus_Giacomo.R")
} 



  
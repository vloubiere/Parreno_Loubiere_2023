setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

if(F)
{
  #-------------------------------------------------------------------------#
  # Available data
  #-------------------------------------------------------------------------#
  #### Process public data and identify/classify ED REs ####
  file.edit("git_epiCancer/vloubiere/subscripts/ED_REs_characterization.R")
  file.edit("git_epiCancer/vloubiere/subscripts/count_promoter_motifs.R")
  
  #-------------------------------------------------------------------------#
  # RNA analyses
  #-------------------------------------------------------------------------#
  #### Process RNA-Seq data ####
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/GFP_counting_RNASeq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/alignment_statistics.R")
  
  #### Bulk analyses ####
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_FPKM_transriptomes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_log2FC_transriptomes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_tissues_vs_transplants.R")
  file.edit("git_epiCancer/vloubiere/subscripts/MA_plots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/alluvial_plot_affected_genes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/upset_plot_overlapping_genes_transriptomes.R")
  
  #### Clustering ####
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_transcriptomes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/figure_transcriptomes_clusters.R")
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_clusters_FC_diff.R")
  
  #### Tables for AMM ####
  file.edit("git_epiCancer/vloubiere/subscripts/tables_AMM.R")
  file.edit("git_epiCancer/vloubiere/subscripts/collect_figures_pdfs_RNA.R")
  
  #-------------------------------------------------------------------------#
  # Cut N run
  #-------------------------------------------------------------------------#
  # Processing and QC
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_QC.R")
  file.edit("git_epiCancer/vloubiere/subscripts/compare_vl_gonza_cutnrun_tracks.R") # Gonzalo's files removed
  
  # Analyses
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_MA_plots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_target_genes_analysis.R")
  
  # To be cleaned
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_clustering_DESeq2.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_clustering_peaks.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_K27_changes_volcano_plot.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_differential_peak_calling.R")
  
  # New runs
  file.edit("git_epiCancer/vloubiere/subscripts/new_cutnruns_PCC_peakCalling.R")
  
  #-------------------------------------------------------------------------#
  # Mutations
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/genomic_SNPs_mutations.R")
  file.edit("git_epiCancer/vloubiere/subscripts/genomic_InDels_mutations.R")
  
  
  #-------------------------------------------------------------------------#
  # MISCELLANEOUS / Casual
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/flowJo.R")
  file.edit("git_epiCancer/vloubiere/subscripts/dac_locus_Giacomo.R")
} 



  
setwd("D:/_R_data/projects/epigenetic_cancer/")

if(F)
{
  #-----------------------------------#
  # METADATA
  #-----------------------------------#
  # Can be used to collect some infos, but the real metadata can be found at "Rdata/raw_metadata_final.txt"
  file.edit("git_epiCancer/vloubiere/subscripts/raw_metadata.R")
  
  #-----------------------------------#
  # Download available data
  #-----------------------------------#
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/download_modEncode_RNASeq.R") # Not working
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/download_Paro_data.R")
  
  #-----------------------------------#
  # Process data
  #-----------------------------------#
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/transcriptomes_processing.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/alignment_statistics.R") #OK
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/MA_plots.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R") #OK
  
  #-----------------------------------#
  # Compute available data
  #-----------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_ATAC_seq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_K9me3_seq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/REs_gene_assignment.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Compute_available_data.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Compute_motifs.R")
  file.edit("git_epiCancer/vloubiere/subscripts/tissue_specificity_score.R")
  
  #-----------------------------------#
  # Kmeans clustering transcriptomes
  #-----------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/clustering.R")
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_heatmaps.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Targeted_heatmaps.R")  
  file.edit("git_epiCancer/vloubiere/subscripts/Gene_ontologies.R")
  file.edit("git_epiCancer/vloubiere/subscripts/String_networks.R")
  
  #-----------------------------------#
  # Further analyses
  #-----------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/motif_wt_activity.R") # Aim was to classify motifs as rep/act. Not very conclusive
  file.edit("git_epiCancer/vloubiere/subscripts/compare_Paro_data.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_transriptomes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_all_transriptomes_FPKMs.R") # In this version, I used FPKMs instead of log2FC
  
  file.edit("git_epiCancer/vloubiere/subscripts/transplantation_FACS_comparison.R")
  
  file.edit("git_epiCancer/vloubiere/subscripts/dac_locus_Giacomo.R")
  
  # file.edit("git_epiCancer/vloubiere/subscripts/transplantation_experiments.R") # Attempt to compute raw FCs for single replicates T5, T8
  
  #--------------------------------------------------------------------------------------------------------------------------------------------#
  # gDNA analyses
  #--------------------------------------------------------------------------------------------------------------------------------------------#
  #-----------------------------------#
  # Processing novogene files
  #-----------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/assemble_VCF_function_files.R")
  file.edit("git_epiCancer/vloubiere/subscripts/melt_data_compute_loh.R")
  file.edit("git_epiCancer/vloubiere/subscripts/fisher_test_alleles.R") # Not working yet
  
  file.edit("git_epiCancer/vloubiere/subscripts/intersection_upsetplots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/pie_chart_mutation_type.R")
  file.edit("git_epiCancer/vloubiere/subscripts/mutations_impact_transcription.R")
  
  #-----------------------------------#
  # Restart from sratch
  #-----------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/download_tano_lgl.R")
  file.edit("git_epiCancer/vloubiere/subscripts/SNPs_calling.R")
  file.edit("git_epiCancer/vloubiere/subscripts/SNP_analysis.R")
  
  file.edit("git_epiCancer/vloubiere/subscripts/polarity_GO_JP.R")
  file.edit("git_epiCancer/vloubiere/subscripts/motifs_enrichment_AMaria.R")
} 



  
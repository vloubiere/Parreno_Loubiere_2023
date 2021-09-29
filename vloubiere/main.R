setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

if(F)
{
  #-------------------------------------------------------------------------#
  # Available data
  #-------------------------------------------------------------------------#
  #### Process public data and identify/classify ED REs ####
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_available_data.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/ATAC_Seq_peak_calling.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/REs_assignment_and_features.R") #OK
  
  #-------------------------------------------------------------------------#
  # RNA analyses
  #-------------------------------------------------------------------------#
  #### Process RNA-Seq data ####
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_processing.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/alignment_statistics.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/generate_final_FC_table.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/generate_tables_for_AMM.R")
  
  #### Bulk analyses ####
  file.edit("git_epiCancer/vloubiere/subscripts/MA_plots.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/alluvial_plot_affected_genes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_transriptomes.R") #OK
  
  #### Clustering DOSE ####
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_dose_transcriptomes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/GO_dose_transcriptomes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/STRING_dose_transcriptomes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/features_dose_transcriptomes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/diff_18_29_dose.R") #OK
  
  #### Clustering epiCancenr ####
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_epiCancer_transcriptomes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/GO_dose_transcriptomes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/STRING_dose_transcriptomes.R") #OK
  
  # NEXT: Check apoptosis pathway in dose timecourse
  
  #### Clustering TS ####
  # NEXT: Cluster TS experiment and check whether simpler groups from the ppt make more sense
  # Work with both clusters and simpler meta clusters??
  # Afeter clustering -> GO, STRING and selected gene heatmaps (make beautiful figure, see heatmaps_selected_genes.R)
  
  
  #### Overlapping genes analyses ####
  file.edit("git_epiCancer/vloubiere/subscripts/upsetPlot_overlapping_genes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/vennDiag_overlapping_genes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/GO_overlapping_genes.R") #OK
  
  #### Kmeans clustering transcriptomes ####
  file.edit("git_epiCancer/vloubiere/subscripts/clustering.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/compute_clusters_enrichment.R") #OK
  
  
  file.edit("git_epiCancer/vloubiere/subscripts/heatmaps_selected_genes.R") #OK
  file.edit("git_epiCancer/vloubiere/subscripts/Targeted_heatmaps.R")  
  file.edit("git_epiCancer/vloubiere/subscripts/Gene_ontologies.R")
  file.edit("git_epiCancer/vloubiere/subscripts/String_networks.R")
  
  #### Further analyses ####
  file.edit("git_epiCancer/vloubiere/subscripts/motif_wt_activity.R") # Aim was to classify motifs as rep/act. Not very conclusive
  file.edit("git_epiCancer/vloubiere/subscripts/compare_Paro_data.R")
  file.edit("git_epiCancer/vloubiere/subscripts/transplantation_FACS_comparison.R")
  
  #-------------------------------------------------------------------------#
  # gDNA analyses
  #-------------------------------------------------------------------------#
  #### Somatic mutations Novogene ####
  file.edit("git_epiCancer/vloubiere/subscripts/somatic_mutations_Novogene.R")
  
  
  
  #### Processing novogene files ####
  file.edit("git_epiCancer/vloubiere/subscripts/assemble_VCF_function_files.R")
  file.edit("git_epiCancer/vloubiere/subscripts/melt_data_compute_loh.R")
  file.edit("git_epiCancer/vloubiere/subscripts/fisher_test_alleles.R") # Not working yet
  file.edit("git_epiCancer/vloubiere/subscripts/intersection_upsetplots.R")
  file.edit("git_epiCancer/vloubiere/subscripts/pie_chart_mutation_type.R")
  file.edit("git_epiCancer/vloubiere/subscripts/mutations_impact_transcription.R")
  
  #### Restart gDNA analyses from scratch ####
  file.edit("git_epiCancer/vloubiere/subscripts/SNPs_calling.R")
  file.edit("git_epiCancer/vloubiere/subscripts/SNP_analysis.R")
  file.edit("git_epiCancer/vloubiere/subscripts/polarity_GO_JP.R")
  file.edit("git_epiCancer/vloubiere/subscripts/motifs_enrichment_AMaria.R")
  
  #-------------------------------------------------------------------------#
  # MISCELLANEOUS / Casual
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/dac_locus_Giacomo.R")
  
  #-------------------------------------------------------------------------#
  # Old/Unused scripts
  #-------------------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/raw_metadata.R") # Can be used to collect some infos, but the real metadata can be found at "Rdata/raw_metadata_final.txt
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/download_Paro_data.R") # Not used anymnore (see raw metadata sheet)
  file.edit("D:/_R_data/projects/epigenetic_cancer/git_epiCancer/vloubiere/subscripts/download_modEncode_RNASeq.R") # Not working
  file.edit("git_epiCancer/vloubiere/subscripts/download_tano_lgl.R") # Not used anymnore (see raw metadata sheet)
  #-------------#
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_all_transriptomes_FPKMs.R") # In this version, I used FPKMs instead of log2FC
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_ATAC_seq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Processing_K9me3_seq.R")
  file.edit("git_epiCancer/vloubiere/subscripts/REs_gene_assignment.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Compute_available_data.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Compute_motifs.R")
  file.edit("git_epiCancer/vloubiere/subscripts/tissue_specificity_score.R")
  #-------------#
  file.edit("git_epiCancer/vloubiere/subscripts/transplantation_experiments.R") # Attempt to compute raw FCs for single replicates T5, T8
} 



  
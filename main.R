setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

if(F)
{
  ###########################################################################
  # Available data
  ###########################################################################
  # Process public data and identify/classify ED REs -----------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/ED_REs_characterization.R")
  file.edit("git_epiCancer/vloubiere/subscripts/count_promoter_motifs.R")
  
  ###########################################################################
  # RNA analyses
  ###########################################################################
  # Process RNA-Seq data ---------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/GFP_counting_RNASeq.R")
  
  # Bulk analyses ----------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/alignment_statistics.R")
  file.edit("git_epiCancer/vloubiere/subscripts/GFP_reads_barplot.R") # Validate genetic system
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_correlations.R") # not used
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_counts_transriptomes.R") # Compare replicates
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_log2FC_transriptomes.R") # Compare conditions
  file.edit("git_epiCancer/vloubiere/subscripts/PCA_tissues_vs_transplants.R") # not used
  file.edit("git_epiCancer/vloubiere/subscripts/MA_plots.R") # Check up down genes
  file.edit("git_epiCancer/vloubiere/subscripts/alluvial_plot_affected_genes.R") # Transitions timecourse
  file.edit("git_epiCancer/vloubiere/subscripts/upset_plot_overlapping_genes_transriptomes.R") # Overlaps between conditions
  file.edit("git_epiCancer/vloubiere/subscripts/heatmap_FC_GOF.R") # heatmap FC genes of interest
  file.edit("git_epiCancer/vloubiere/subscripts/GO_up_down_genes.R") # For each RNAi cdition, up/down PRC1+/1 GOs
  
  ###########################################################################
  # Cut N run
  ###########################################################################
  # Processing and QC ------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_processing.R")
  file.edit("git_epiCancer/vloubiere/subscripts/ecdysone_cutnrun_processing.R")
  
  # QC checks --------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/spike-in_perc.R") # Not used
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_screenshots.R") # Just to check
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_PCC_replicates.R") # PCC reps
  
  # Analyses ---------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_MA_plots.R")
  
  # Not used ---------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/compare_vl_gonza_cutnrun_tracks.R") # Gonzalo's files removed
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_target_genes_analysis.R") # Not finished
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_clustering_DESeq2.R") # Not used
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_clustering_peaks.R") # Not used
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_K27_changes_volcano_plot.R") # Not used
  file.edit("git_epiCancer/vloubiere/subscripts/cutnrun_differential_peak_calling.R") # Not used
  
  ###########################################################################
  # Clustering
  ###########################################################################
  # Clustering -------------------------------------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/clustering_transcriptomes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/figure_transcriptomes_clusters.R")
  
  # Quantif features / cluster (PRC1+/-) -----------------------------------#
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_clusters_FC_diff.R")
  file.edit("git_epiCancer/vloubiere/subscripts/transcriptomes_clusters_FC_development.R") # Target genes are dev genes? Not convincing
  
  ###########################################################################
  # Reverting vs non reverting 
  ###########################################################################
  file.edit("git_epiCancer/vloubiere/subscripts/define_RECOVERY_&_NORECOVERY_genes.R")
  file.edit("git_epiCancer/vloubiere/subscripts/compare_RECOVERY_vs_transcriptome_clusters_2&5.R")
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_screenshots.R")
  # Promoters average tracks
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_average_tracks_PH18.R")
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_average_tracks_compare_conditions.R")
  # Motifs proximal REs
  file.edit("git_epiCancer/vloubiere/subscripts/assign_ATAC_Seq_peaks_to_genes.R") # Also includes promoters and PH peaks
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_REs_motifs.R")
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_REs_average_tracks_PH18.R")
  # FC during development (not used)
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_FC_development.R")
  file.edit("git_epiCancer/vloubiere/subscripts/RECOVERY_vs_NOT_FPKM_FLYBASE_development.R")
  file.edit("git_epiCancer/vloubiere/subscripts/Overlap_with_APF_induced_genes.R")
  
  ###########################################################################
  # Cluster 2 vs cluster 5
  ###########################################################################
  file.edit("git_epiCancer/vloubiere/subscripts/define_cl2_&_cl5_target_genes.R")
  # Promoters average tracks
  file.edit("git_epiCancer/vloubiere/subscripts/cl2_vs_cl5_average_tracks_PH18.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cl2_vs_cl5_average_tracks_compare_conditions.R")
  # Motifs proximal REs (in this case, PH peaks!)
  file.edit("git_epiCancer/vloubiere/subscripts/cl2_vs_cl5_REs_motifs.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cl2_vs_cl5_REs_average_tracks_PH18.R")
  file.edit("git_epiCancer/vloubiere/subscripts/cl2_&_cl5_target_genes_ph_peaks_density.R")
  
  ###########################################################################
  # Tables / Figures for AMM  ----------------------------------------------#
  ###########################################################################
  file.edit("git_epiCancer/vloubiere/subscripts/tables_AMM.R")
  file.edit("git_epiCancer/vloubiere/subscripts/update_dropbox.R")
  file.edit("git_epiCancer/vloubiere/subscripts/update_dropbox.R")
  
  ###########################################################################
  # Mutations
  ###########################################################################
  file.edit("git_epiCancer/vloubiere/subscripts/genomic_SNPs_mutations.R")
  file.edit("git_epiCancer/presentation.Rhtml")
  
  ###########################################################################
  # MISCELLANEOUS / Casual
  ###########################################################################
  file.edit("git_epiCancer/vloubiere/subscripts/flowJo.R")
  file.edit("git_epiCancer/vloubiere/subscripts/dac_locus_Giacomo.R")
}
  
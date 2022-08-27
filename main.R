setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

if(F)
{
  ###########################################################################
  # Available data
  ###########################################################################
  # Process public data and identify/classify ED REs -----------------------#
  file.edit("git_epiCancer/subscripts/ED_REs_characterization.R")
  file.edit("git_epiCancer/subscripts/count_promoter_motifs.R")
  
  ###########################################################################
  # RNA analyses
  ###########################################################################
  # Process RNA-Seq data ---------------------------------------------------#
  file.edit("git_epiCancer/subscripts/benchmark_transcriptomes_processing.R")
  file.edit("git_epiCancer/subscripts/transcriptomes_processing.R")
  file.edit("git_epiCancer/subscripts/GFP_counting_RNASeq.R")
  
  # Bulk analyses ----------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/alignment_statistics.R")
  file.edit("git_epiCancer/subscripts/GFP_reads_barplot.R") # Validate genetic system
  file.edit("git_epiCancer/subscripts/transcriptomes_correlations.R") # not used
  file.edit("git_epiCancer/subscripts/PCA_counts_transriptomes.R") # Compare replicates
  file.edit("git_epiCancer/subscripts/PCA_log2FC_transriptomes.R") # Compare conditions
  file.edit("git_epiCancer/subscripts/PCA_tissues_vs_transplants.R") # not used
  file.edit("git_epiCancer/subscripts/MA_plots.R") # Check up down genes
  file.edit("git_epiCancer/subscripts/alluvial_plot_affected_genes.R") # Transitions timecourse
  file.edit("git_epiCancer/subscripts/upset_plot_overlapping_genes_transriptomes.R") # Overlaps between conditions
  file.edit("git_epiCancer/subscripts/heatmap_FC_GOF.R") # heatmap FC genes of interest
  file.edit("git_epiCancer/subscripts/GO_up_down_genes.R") # For each RNAi cdition, up/down PRC1+/1 GOs
  
  ###########################################################################
  # Cut N run
  ###########################################################################
  # Processing and QC ------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/cutnrun_processing.R")
  file.edit("git_epiCancer/subscripts/ecdysone_cutnrun_processing.R")
  
  # QC checks --------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/spike-in_perc.R") # Not used
  file.edit("git_epiCancer/subscripts/cutnrun_screenshots.R") # Just to check
  file.edit("git_epiCancer/subscripts/cutnrun_PCC_replicates.R") # PCC reps
  
  # Analyses ---------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/cutnrun_MA_plots.R")
  
  # Not used ---------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/compare_vl_gonza_cutnrun_tracks.R") # Gonzalo's files removed
  file.edit("git_epiCancer/subscripts/cutnrun_target_genes_analysis.R") # Not finished
  file.edit("git_epiCancer/subscripts/cutnrun_clustering_DESeq2.R") # Not used
  file.edit("git_epiCancer/subscripts/cutnrun_clustering_peaks.R") # Not used
  file.edit("git_epiCancer/subscripts/cutnrun_K27_changes_volcano_plot.R") # Not used
  file.edit("git_epiCancer/subscripts/cutnrun_differential_peak_calling.R") # Not used
  
  ###########################################################################
  # Clustering
  ###########################################################################
  # Clustering -------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/clustering_transcriptomes.R") # SOM
  file.edit("git_epiCancer/subscripts/make_GENEs_data_table.R") # For each gene, retrieve RNA-Seq and compute features 
  file.edit("git_epiCancer/subscripts/make_REs_data_table.R") # Assign ATAC, TSS and PH peaks to closest prom -> motif +/-250bp
  
  # Clusters analysis
  file.edit("git_epiCancer/subscripts/compute_clusters_extra_features.R") # Network, GO, motifs and 
  file.edit("git_epiCancer/subscripts/figure_transcriptomes_clusters.R") # Clusters Fig. with motifs and GO analysis
  
  # Quantif features / cluster (PRC1+/-) -----------------------------------#
  file.edit("git_epiCancer/subscripts/transcriptomes_clusters_FC_diff.R") # RNA-Seq/CHIP signal per cluster +/- PRC1
  file.edit("git_epiCancer/subscripts/transcriptomes_clusters_FC_development.R") # Target genes are dev genes? Not convincing
  
  ###########################################################################
  # Cluster 2 vs cluster 5
  ###########################################################################
  # Promoters average tracks
  file.edit("git_epiCancer/subscripts/cl2_vs_cl5_average_tracks_PH18.R")
  file.edit("git_epiCancer/subscripts/cl2_vs_cl5_average_tracks_compare_conditions.R")
  # Motifs proximal REs (in this case, PH peaks!)
  file.edit("git_epiCancer/subscripts/cl2_vs_cl5_REs_motifs.R")
  file.edit("git_epiCancer/subscripts/cl2_vs_cl5_REs_average_tracks_PH18.R")
  file.edit("git_epiCancer/subscripts/cl2_&_cl5_target_genes_ph_peaks_density.R")
  
  ###########################################################################
  # Reverting vs non reverting 
  ###########################################################################
  file.edit("git_epiCancer/subscripts/define_RECOVERY_&_NORECOVERY_genes.R")
  file.edit("git_epiCancer/subscripts/compare_RECOVERY_vs_transcriptome_clusters_2&5.R")
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_screenshots.R")
  # Promoters average tracks
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_average_tracks_PH18.R")
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_average_tracks_compare_conditions.R")
  # Motifs proximal REs
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_REs_motifs.R")
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_REs_average_tracks_PH18.R")
  # FC during development (not used)
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_FC_development.R")
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_FPKM_FLYBASE_development.R")
  file.edit("git_epiCancer/subscripts/Overlap_with_APF_induced_genes.R")
  
  ###########################################################################
  # Tables / Figures for AMM  ----------------------------------------------#
  ###########################################################################
  file.edit("git_epiCancer/subscripts/tables_AMM.R")
  file.edit("git_epiCancer/subscripts/update_dropbox.R")
  file.edit("git_epiCancer/epiCancer_presentation.Rmd")
  
  ###########################################################################
  # Mutations
  ###########################################################################
  file.edit("git_epiCancer/subscripts/genomic_SNPs_mutations.R")
  
  ###########################################################################
  # MISCELLANEOUS / Casual
  ###########################################################################
  file.edit("git_epiCancer/subscripts/flowJo.R")
  file.edit("git_epiCancer/subscripts/dac_locus_Giacomo.R")
}
  
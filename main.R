setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

if(F)
{
  ###########################################################################
  # DATA processing
  ###########################################################################
  # Trim/align/count/diff. analysis. ---------------------------------------#
  file.edit("git_epiCancer/subscripts/transcriptomes_processing.R")
  file.edit("git_epiCancer/subscripts/cutnrun_processing.R")
  file.edit("git_epiCancer/subscripts/ATAC_Seq_processing.R")
  file.edit("git_epiCancer/subscripts/ecdysone_cutnrun_processing.R")
  file.edit("git_epiCancer/subscripts/make_REs_data_table.R") # Assign TSS/ATAC/PH peaks to closest prom -> motif +/-250bp
  
  ###########################################################################
  # RNA analyses
  ###########################################################################
  file.edit("git_epiCancer/subscripts/Alignment_statistics_RNA.R") #QC
  file.edit("git_epiCancer/subscripts/GFP_reads_barplot.R") # Validate genetic system
  file.edit("git_epiCancer/subscripts/transcriptomes_correlations.R") # PCC between replicates
  file.edit("git_epiCancer/subscripts/PCA_counts_transriptomes.R") # Compare replicates
  file.edit("git_epiCancer/subscripts/PCA_log2FC_transriptomes.R") # Compare conditions
  file.edit("git_epiCancer/subscripts/MA_plots.R") # Check up down genes
  file.edit("git_epiCancer/subscripts/alluvial_plot_affected_genes.R") # Transitions timecourse
  file.edit("git_epiCancer/subscripts/upset_plot_overlapping_genes_transriptomes.R") # Overlaps between conditions
  file.edit("git_epiCancer/subscripts/compare_GFP_noGFP.R") # Overlaps between systems
  file.edit("git_epiCancer/subscripts/heatmap_FC_GOF.R") # heatmap FC genes of interest
  file.edit("git_epiCancer/subscripts/GO_up_down_genes.R") # For each RNAi cdition, up/down PRC1+/1 GOs
  
  file.edit("git_epiCancer/subscripts/PCA_tissues_vs_transplants.R") # not used
  
  ###########################################################################
  # Cut N run
  ###########################################################################
  file.edit("git_epiCancer/subscripts/cutnrun_PCC_replicates.R") # PCC reps
  file.edit("git_epiCancer/subscripts/cutnrun_screenshots.R") # Just to check
  file.edit("git_epiCancer/subscripts/Alignment_statistics_cutnrun.R") # Not used
  file.edit("git_epiCancer/subscripts/spike-in_perc.R") # Not used
  
  # Analyses ---------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/cutnrun_MA_plots.R") # based on peaks
  file.edit("git_epiCancer/subscripts/cutnrun_interesect_diff_regions.R") # Based on collapsed diff peaks -> overlap?
  
  ###########################################################################
  # Clustering
  ###########################################################################
  file.edit("git_epiCancer/subscripts/clustering_transcriptomes.R") # SOM
  file.edit("git_epiCancer/subscripts/make_GENEs_data_table.R") # For each gene, retrieve RNA-Seq and compute features 
  file.edit("git_epiCancer/subscripts/compute_clusters_extra_features.R") # Network, GO, motifs and 
  file.edit("git_epiCancer/subscripts/figure_transcriptomes_clusters.R") # Clusters Fig. with motifs and GO analysis
  file.edit("git_epiCancer/subscripts/network_per_clusters.R") # Clusters networks
  file.edit("git_epiCancer/subscripts/transcriptomes_clusters_FC_diff.R") # RNA-Seq/CHIP signal per cluster +/- PRC1
  file.edit("git_epiCancer/subscripts/cluster_genes_tissue_specificity.R") # Based on encode method
  
  ###########################################################################
  # Reverting vs non reverting 
  ###########################################################################
  file.edit("git_epiCancer/subscripts/Characterization_recovery_groups.R") # Overlap RNA clusters and FC
  file.edit("git_epiCancer/subscripts/Recovery_groups_extra_features.R") # GO, network and motifs recov vs no recov
  file.edit("git_epiCancer/subscripts/recovery_average_tracks.R") # average tracks in control and between cditions
  file.edit("git_epiCancer/subscripts/recovery_ChIP_FC_cditions.R") # ChIP FC per condition
  file.edit("git_epiCancer/subscripts/recovery_ChIP_percentage_cditions.R") # Similar to previous one, but with % of WT coverage
  file.edit("git_epiCancer/subscripts/recovery_FPKM_cditions.R")
  file.edit("git_epiCancer/subscripts/recovery_K36me3_enrich.R")
  
  # Not used --------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_screenshots.R")
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_FC_development.R")
  file.edit("git_epiCancer/subscripts/Overlap_with_APF_induced_genes.R")
  
  ###########################################################################
  # Genomic DNA
  ###########################################################################
  file.edit("git_epiCancer/subscripts/gDNA_final_table.R")
  file.edit("git_epiCancer/subscripts/gDNA_overlaps_conditions.R")
  file.edit("git_epiCancer/subscripts/gDNA_mut_counts.R")
  file.edit("git_epiCancer/subscripts/gDNA_allele_frequency.R")
  file.edit("git_epiCancer/subscripts/gDNA_hotspots.R")
  
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
  
# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

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
  
  ###########################################################################
  # Clustering and data tables
  ###########################################################################
  file.edit("git_epiCancer/subscripts/clustering_transcriptomes.R") # SOM
  file.edit("git_epiCancer/subscripts/compute_clusters_extra_features.R") # Network, GO, motifs and STRING network
  file.edit("git_epiCancer/subscripts/make_GENEs_data_table.R") # For each gene, retrieve RNA-Seq and compute features 
  file.edit("git_epiCancer/subscripts/gDNA_final_table.R")
  file.edit("git_epiCancer/subscripts/make_REs_data_table.R") # not used
  
  ###########################################################################
  # Main figures
  ###########################################################################
  # Figure 1
  file.edit("git_epiCancer/subscripts/barplot_WB_quantif.R") # Not used
  file.edit("git_epiCancer/subscripts/violin_plot_phKD_sizes.R")
  file.edit("git_epiCancer/subscripts/mutations_features.R")
  file.edit("git_epiCancer/subscripts/exonic_mutations_features.R")
  file.edit("git_epiCancer/subscripts/mutated_genes_overlaps.R")
  
  # Supplementary
  file.edit("git_epiCancer/subscripts/gDNA_GO_mutated_genes.R") # No GO enriched for genes mutated in at least two conditions
  
  # Figure 2
  file.edit("git_epiCancer/subscripts/MA_plots.R") # Check up down genes
  file.edit("git_epiCancer/subscripts/alluvial_plot_affected_genes.R") # Transitions timecourse
  file.edit("git_epiCancer/subscripts/figure_transcriptomes_clusters.R") # Clusters Fig. with network and motifs
  file.edit("git_epiCancer/subscripts/clusters_GOs.R") # Clusters GO analysis
  file.edit("git_epiCancer/subscripts/heatmap_FC_JAKSTAT_JNK.R") # heatmap FC genes of interest
  file.edit("git_epiCancer/subscripts/violin_plot_zfh1KD_sizes.R")
  
  # Figure 3
  file.edit("git_epiCancer/subscripts/boxplots_FPKM_CUTNRUN_CHIP.R") # Quantif signal
  file.edit("git_epiCancer/subscripts/barplot_PH.R") # Count overlap with peaks
  file.edit("git_epiCancer/subscripts/screenshots.R")
  file.edit("git_epiCancer/subscripts/MA_plots_K27.R") # based on peaks
  
  # Figure 4
  file.edit("git_epiCancer/subscripts/boxplot_FPKMs_PHbound_genes.R") # average tracks in control and between cditions
  file.edit("git_epiCancer/subscripts/recovery_average_tracks.R") # average tracks in control and between cditions
  file.edit("git_epiCancer/subscripts/recovery_average_tracks_PH_peaks.R") # average tracks in control and between cditions
  file.edit("git_epiCancer/subscripts/recovery_motifs_enrichment_promoters.R") # Motifs enrichment promoters
  file.edit("git_epiCancer/subscripts/recovery_motifs_enrichment_PH_peaks.R") # Motifs enrichment PH peaks
  file.edit("git_epiCancer/subscripts/recovery_motifs_icisTarget_ATAC_peaks.R") # Motifs enrichment at ATAC-Seq peaks
  file.edit("git_epiCancer/subscripts/recovery_motifs_icisTarget_genes_list.R") # Not used -> Motifs enrichment at genes
  file.edit("git_epiCancer/subscripts/recovery_motifs_icisTarget_PH_peaks.R") # Not used -> Motifs enrichment at PH peaks
  
  
  # Extended data 2
  file.edit("git_epiCancer/subscripts/violin_plot_PscKD_sizes.R")
  
  # Extended data 5
  file.edit("git_epiCancer/subscripts/zfh1_screenshot.R")
  file.edit("git_epiCancer/subscripts/percentage_PRC1_bound_genes.R") # Check percentage of PRC1 bound genes per cluster
  file.edit("git_epiCancer/subscripts/transcriptomes_clusters_FC_diff.R") # RNA-Seq/CHIP signal per cluster +/- PRC1
  file.edit("git_epiCancer/subscripts/MA_plot_K118Ub.R") # based on peaks
  file.edit("git_epiCancer/subscripts/recovery_average_tracks_K27_K118.R") # average tracks in control and between cditions
  # For PH quantif -> see figure 3
  
  # Extended data 6
  file.edit("git_epiCancer/subscripts/GFP_reads_barplot.R") # Validate genetic system
  
  ###########################################################################
  # RNA analyses
  ###########################################################################
  file.edit("git_epiCancer/subscripts/Alignment_statistics_RNA.R") #QC
  file.edit("git_epiCancer/subscripts/transcriptomes_correlations.R") # PCC between replicates
  file.edit("git_epiCancer/subscripts/PCA_counts_transriptomes.R") # Compare replicates
  file.edit("git_epiCancer/subscripts/PCA_log2FC_transriptomes.R") # Compare conditions
  file.edit("git_epiCancer/subscripts/chisq_compare_GFP_noGFP_up_down_genes.R") # Overlaps between systems
  file.edit("git_epiCancer/subscripts/boxplot_FC_up_down_genes.R") # Check up down genes
  file.edit("git_epiCancer/subscripts/upset_plot_overlapping_genes_transriptomes.R") # Overlaps between conditions
  file.edit("git_epiCancer/subscripts/heatmap_FC_GOF.R") # heatmap FC genes of interest
  file.edit("git_epiCancer/subscripts/GO_up_down_genes.R") # For each RNAi cdition, up/down PRC1+/1 GOs
  
  file.edit("git_epiCancer/subscripts/PCA_tissues_vs_transplants.R") # not used
  
  ###########################################################################
  # RNA clusters analyses
  ###########################################################################
  file.edit("git_epiCancer/subscripts/network_per_clusters.R") # Clusters networks
  file.edit("git_epiCancer/subscripts/cluster_genes_tissue_specificity.R") # Based on encode method
  
  ###########################################################################
  # Cut N run
  ###########################################################################
  file.edit("git_epiCancer/subscripts/cutnrun_PCC_replicates.R") # PCC reps
  file.edit("git_epiCancer/subscripts/Alignment_statistics_cutnrun.R") # Read counts
  
  # Analyses ---------------------------------------------------------------#
  file.edit("git_epiCancer/subscripts/cutnrun_MA_plots.R") # based on peaks
  file.edit("git_epiCancer/subscripts/cutnrun_interesect_diff_regions.R") # Based on collapsed diff peaks -> overlap?
  
  file.edit("git_epiCancer/subscripts/cutnrun_screenshots.R") # Not used
  file.edit("git_epiCancer/subscripts/spike-in_perc.R") # Not used
  
  ###########################################################################
  # Reverting vs non reverting 
  ###########################################################################
  file.edit("git_epiCancer/subscripts/recovery_FC_conditions.R") # FC between conditions
  file.edit("git_epiCancer/subscripts/recovery_overlap_RNA_clusters.R") # Overlap RNA clusters
  file.edit("git_epiCancer/subscripts/recovery_chromhmm_class.R") # Chromhmm classes
  file.edit("git_epiCancer/subscripts/recovery_overlap_K27_K118_domains.R") # Overlap with K27 and K118 domains
  file.edit("git_epiCancer/subscripts/recovery_size_K27_domains.R") # Size of ovarlapping K27 domains
  file.edit("git_epiCancer/subscripts/recovery_FPKM_cditions.R") # Check FPKMs recovery genes
  file.edit("git_epiCancer/subscripts/recovery_GO.R") # GO
  file.edit("git_epiCancer/subscripts/recovery_STRING_network.R") # Network
  file.edit("git_epiCancer/subscripts/recovery_motifs_enrichment.R") # Motifs enrichment
  file.edit("git_epiCancer/subscripts/recovery_motifs_counts.R") # Motifs coutns histograms
  file.edit("git_epiCancer/subscripts/recovery_motifs_specific_proteins.R") # Heatmap ZEB STAT motifs
  file.edit("git_epiCancer/subscripts/recovery_average_tracks_neoPRC1.R") # average tracks in control and between cditions
  file.edit("git_epiCancer/subscripts/recovery_ChIP_FC_cditions.R") # ChIP FC per condition
  file.edit("git_epiCancer/subscripts/recovery_K36me3_enrich.R") # No recovery genes abortive transcription?
  file.edit("git_epiCancer/subscripts/recovery_K27Ac_enrich.R") # Recovery genes have particular chromatin?
  file.edit("git_epiCancer/subscripts/recovery_PH_levels_vs_FC.R") # PH levels vs FC
  file.edit("git_epiCancer/subscripts/recovery_PH_K27me3_peaks_overlaps_between_cditions.R") # PH levels vs FC
  
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_screenshots.R") # Not used
  file.edit("git_epiCancer/subscripts/RECOVERY_vs_NOT_FC_development.R") # Not used
  file.edit("git_epiCancer/subscripts/Overlap_with_APF_induced_genes.R") # Not used
  
  ###########################################################################
  # Genomic DNA
  ###########################################################################
  file.edit("git_epiCancer/subscripts/gDNA_overlaps_conditions.R")
  file.edit("git_epiCancer/subscripts/gDNA_mut_counts.R")
  file.edit("git_epiCancer/subscripts/gDNA_allele_frequency.R")
  file.edit("git_epiCancer/subscripts/gDNA_overlaps_mutations_shared_between_several_cditions.R")
  file.edit("git_epiCancer/subscripts/gDNA_mut_features.R")
  file.edit("git_epiCancer/subscripts/gDNA_exonic_mut_features.R")
  file.edit("git_epiCancer/subscripts/gDNA_genes_overlaps.R")
  file.edit("git_epiCancer/subscripts/gDNA_table_mutated_genes.R")
  file.edit("git_epiCancer/subscripts/gDNA_mutations_impact_FC.R")
  file.edit("git_epiCancer/subscripts/gDNA_hotspots.R")
  
  ###########################################################################
  # Figures ----------------------------------------------------------------#
  ###########################################################################
  file.edit("git_epiCancer/epiCancer_presentation.Rmd")
  file.edit("git_epiCancer/subscripts/update_dropbox.R")
  
  
  
  
  ###########################################################################
  # MISCELLANEOUS / Casual
  ###########################################################################
  file.edit("git_epiCancer/subscripts/flowJo.R")
  file.edit("git_epiCancer/subscripts/dac_locus_Giacomo.R")
}
  
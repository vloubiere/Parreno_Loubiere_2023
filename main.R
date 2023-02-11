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
}
  
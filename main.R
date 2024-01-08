# Data processing ------------------------------------------------------------------------
# Transcriptomes
file.edit("subscripts/transcriptomes_processing.R")
file.edit("subscripts/review_transcriptomes_processing.R")
# Cut&run/ChIP-Seq
file.edit("subscripts/cutnrun_processing.R")
file.edit("subscripts/review_cutnrun_processing.R")
# ATAC-Seq
file.edit("subscripts/review_ATAC_processing.R")

# Make data tables ----------------------------------------------------------------
# gDNA
file.edit("subscripts/review_gDNA_final_table.R")
file.edit("subscripts/review_gDNA_SVs_final_table.R")
file.edit("subscripts/review_gDNA_CNVs_final_table.R")
# cutNrun
file.edit("subscripts/define_PcG_domains.R")
# Transcriptomes
file.edit("subscripts/clustering_transcriptomes.R") # SOM
file.edit("subscripts/make_GENEs_data_table.R")
# ATAC-Seq
file.edit("subscripts/review_ATAC_clustering.R")
file.edit("subscripts/make_ATAC_data_table.R")
file.edit("subscripts/make_rescue_ATAC_data_table.R")

# FIGURES -------------------------------------------------------------------------
# Figure 1 ----
# Disc sizes
file.edit("subscripts/violin_plot_phKD_sizes.R")
# gDNA mutations
file.edit("subscripts/review_allelic_ratio_SNP_InDels.R") # Allele freq SNP
file.edit("subscripts/review_genetic_variants_tumor_fraction.R") # Number of tumors samples
file.edit("subscripts/review_features_SNP_InDels.R") # Features of shared vs control SNP
# DNA repair
file.edit("subscripts/DNA_repair_chiolo.R")

# Figure 2 ----
# Transcriptomes diff analysis
file.edit("subscripts/alluvial_plot_affected_genes.R") # Transitions timecourse
# Gene clusters
file.edit("subscripts/review_transcriptomes_clusters_heatmap.R") # Heatmap clusters
# PcG target genes enrichment
file.edit("subscripts/review_RNA_Clusters_PcG_binding.R") # Barplot over-representation PcG binding
# GO
file.edit("subscripts/review_compute_clusters_GOs.R") # Clusters GO analysis
file.edit("subscripts/review_clusters_GOs.R")
# JAK-STAT pathway
file.edit("subscripts/review_heatmap_FC_JAKSTAT.R") # JAK-STAT pathway

# Figure 3 ----
# Transcription
file.edit("subscripts/review_boxplots_FPKM_CUTNRUN_CHIP.R") # Quantif signal
# PcG binding (peaks/domains)
file.edit("subscripts/review_PcG_binding_irreversible_reversible.R") # Barplots overlap PcG marks rev/irrev
file.edit("subscripts/review_K27Ac_binding_irreversible_reversible.R") # Barplots overlap K27Ac marks rev/irrev
# Screenshot
file.edit("subscripts/review_screenshots.R") # Example irreversible reversible loci
# Classify K27me3 domains
file.edit("subscripts/review_classify_K27_domains_reversible_irreversible.R") # Barplot K27 domains containing rev/irrev
# MA plots
file.edit("subscripts/MA_plots_HTMs.R") # MA plots CUTNRUN peaks FC
# HTMs FC
file.edit("subscripts/review_K27_domains_FC_reversible_irreversible.R")

# Figure 4 ----
# ATAC-Seq clustering
file.edit("subscripts/review_figure_ATAC_clustering.R")
# ATAC-Seq vs RNA clusters
file.edit("subscripts/review_barplot_ATAC_peaks_clusters_association_RNA_clusters.R")
file.edit("subscripts/review_ATAC_peaks_TSS_distance.R")
# Screenshots examples
file.edit("subscripts/review_irreversible_reversible_ATAC_peaks_screenshots.R")
# Motif analysis
file.edit("subscripts/review_ATAC_peaks_iCisTarget_enrichments.R")
file.edit("subscripts/review_predict_ATAC_motifs_LASSO_and_lm.R") # Select relevant motifs and lm
file.edit("subscripts/review_motif_predict_ATAC_FC_lm.R")
file.edit("subscripts/review_ATAC_FC_vs_motif_counts.R")

# Figure 5 ----
# ED sizes
file.edit("subscripts/violin_plot_rescue_RNAi_sizes.R")
# Rescue transcriptomes
file.edit("subscripts/review_alluvial_plot_affected_rescue.R")
# Rescue ATAC-Seq
file.edit("subscripts/review_alluvial_plot_affected_ATAC_rescue.R")
file.edit("subscripts/review_ATAC_FC_rescue_vs_motif_counts.R")
file.edit("subscripts/review_FC_closest_gene_rescue_ATAC.R")
file.edit("subscripts/review_GO_enrichment_rescued_genes.R")

# SUPPLEMENTARY --------------------------------------------------------------------
# Extended data 1 ----
file.edit("subscripts/review_PH_WB.R")
# Extended data 1 ----
file.edit("subscripts/review_RNASeq_cluster_Stat92E_RNAi.R")
# Extended data 2 ----
file.edit("subscripts/violin_plot_PscKD_sizes.R")
# Extended data 3 ----
file.edit("subscripts/review_genetic_variants_overlaps.R")
# Extended data 4 ----
file.edit("subscripts/PCA_counts_transriptomes.R")
file.edit("subscripts/heatmap_genes_of_interest.R")
# Extended data 5 ----
file.edit("subscripts/upset_plot_overlapping_genes_transriptomes.R")
file.edit("subscripts/review_clusters_GOs_full_table.R")
# Extended data 6 ----
file.edit("subscripts/review_average_tracks_promoters_irrev_rev.R")
file.edit("subscripts/overlap_PcG_peaks_and_domains.R")
file.edit("subscripts/review_barplot_classif_K27_domains.R")
file.edit("subscripts/review_average_tracks_PH_peaks_irrev_rev.R")
# Extended data 7 ----
file.edit("subscripts/discs_sizes_transient_rescue.R")
# Extended data 8 ----
file.edit("subscripts/GFP_reads_barplot.R")
file.edit("subscripts/chisq_compare_GFP_noGFP_up_down_genes.R")

# Extended data tables ----
file.edit("subscripts/review_extended_data_table_2.R")
file.edit("subscripts/review_extended_data_table_3.R")
file.edit("subscripts/review_extended_data_table_4.R")
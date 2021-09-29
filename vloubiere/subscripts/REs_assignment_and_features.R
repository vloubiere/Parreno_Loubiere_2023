setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)
require(TFBSTools)

peaks <- fread("db/peaks/ATAC_seq_peaks_ED.txt")

genes <- import("/mnt/d/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
genes <- as.data.table(genes)

#-----------------------------------------#
# Assemble RE table containing TSSs and non-TSS ATAC-Seq peaks
#-----------------------------------------#
TSS <- rbind(genes[type %in% c("mRNA", "ncRNA"), .(seqnames= paste0("chr", seqnames), 
                                                   start= ifelse(strand=="+", start, end),
                                                   end= ifelse(strand=="+", start, end),
                                                   strand, 
                                                   FBgn= gene_id,
                                                   type= "TSS")])
nonTSS <- TSS[peaks, ifelse(.N>0, F, T), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
nonTSS <- peaks[(nonTSS)]
nonTSS$FBgn <- TSS[nonTSS, ifelse(.N>0, 
                                  FBgn[which.min(abs(start-max_coor))],
                                  as.character(NA)), 
                   .EACHI, 
                   on= "seqnames"]$V1
nonTSS[, type:="ENHANCER"]
REs <- rbind(unique(TSS), 
             nonTSS[, .(seqnames, 
                        start= max_coor, 
                        end= max_coor, 
                        strand= "*", 
                        FBgn, 
                        type)])
REs <- na.omit(REs[seqnames %in% c("chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX", "chrY")])

#-----------------------------------------#
# Add gene coordinates
#-----------------------------------------#
REs[genes[type=="gene"], gene_coor:= paste0("chr", 
                                            i.seqnames, 
                                            ":", 
                                            i.start, 
                                            "-", 
                                            i.end,
                                            ":",
                                            i.strand), on= "FBgn==gene_id"]
setcolorder(REs, c("seqnames", "start", "end", "strand", "FBgn", "gene_coor"))

#-----------------------------------------#
# Add ATAC-Seq and PcG binding binary cols
#-----------------------------------------#
REs$open <- peaks[REs, ifelse(.N>0, T, F), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1
PcG <- as.data.table(get(load("external_data/SA2020_cl.list"))$ED_summits)
colnames(PcG)[1:3] <- c("seqnames", "start", "end")
PcG[, start:= start-500]
PcG[, end:= end+500]
REs$PRC1_bound <- PcG[REs, ifelse(.N>0, T, F), .EACHI, on= c("seqnames", "start<=end", "end>=start")]$V1

#-----------------------------------------#
# Add HTMs
#-----------------------------------------#
ext <- REs[, .(seqnames, 
               start= start-500,
               end= end+500)]
REs$ATAC_score <- vl_bw_coverage(ext,
                                 "/mnt/d/_R_data/projects/public_data/dm6/bw/merged_tracks/ATAC_ED_merged.bw")$score
REs$PH_score <- vl_bw_coverage(ext,
                               "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/PH_ED_merge.bw")$score
REs$SUZ12_score <- vl_bw_coverage(ext,
                                  "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/SUZ12_ED_merge.bw")$score
REs$H3K4me1_score <- vl_bw_coverage(ext,
                                    "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/H3K4me1_ED_merge.bw")$score
REs$H3K4me3_score <- vl_bw_coverage(ext,
                                    "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/H3K4me3_ED_merge.bw")$score
REs$H3K27ac_score <- vl_bw_coverage(ext,
                                    "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/H3K27Ac_ED_merge.bw")$score
REs$H3K27me3_score <- vl_bw_coverage(ext,
                                     "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/H3K27me3_ED_merge.bw")$score
REs$H3K27me2_score <- vl_bw_coverage(ext,
                                     "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/H3K27me2_ED_merge.bw")$score
REs$H2AK118Ub_score <- vl_bw_coverage(ext,
                                      "/mnt/d/_R_data/projects/public_data/dm6/bw/SA_2020/H2AK118Ub_ED_merge.bw")$score

#-----------------------------------------#
# add IDS
#-----------------------------------------#
REs[, RE_ID:= .I]
setcolorder(REs, c("seqnames", "start", "end", "strand", "RE_ID"))

#-----------------------------------------#
# Chromatin type clustering
#-----------------------------------------#
mat <- REs[, H3K4me1_score:H2AK118Ub_score]
mat <- apply(mat, 
             2, 
             function(x) 
             {
              lim <- quantile(x, c(0.025, 0.975), na.rm= T)
              x[x<lim[1]] <- lim[1]
              x[x>lim[2]] <- lim[2]
              x <- scale(x)
              return(x)
             })
rownames(mat) <- REs$RE_ID
mat <- na.omit(mat)
boxplot(mat)

pdf("pdf/clustering_REs.pdf")
cl <- vl_heatmap(mat, 
                 kmeans_k = 4, 
                 col = c("royalblue2", "yellow"))
dev.off()

cl[rcl==1, cl_name:= "Polycomb"]
cl[rcl==2, cl_name:= "Null"]
cl[rcl==3, cl_name:= "Enhancer"]
cl[rcl==4, cl_name:= "Tss"]
cl[, row:= as.integer(row)]
REs[cl, HTM_cl:= i.cl_name, on= "RE_ID==row"]
setcolorder(REs, c("seqnames", "start", "end", "strand", "RE_ID", "HTM_cl"))

#-----------------------------------------#
# Motif enrichment
#-----------------------------------------#
sel <- vl_Dmel_motifs_DB$metadata[!is.na(vl_Dmel_motifs_DB$metadata$Dmel) & # Associated to a known TF
                                    vl_Dmel_motifs_DB$metadata$X..motif_collection_name %in% # From a relevant DB
                                    c("flyfactorsurvey", "bergman", "jaspar", "idmmpmm", "cisbp"), "motif_name"]
sel <- which(name(vl_Dmel_motifs_DB$All_pwms_log_odds) %in% sel)
hit <- matchMotifs(vl_Dmel_motifs_DB$All_pwms_log_odds[sel], 
                   resize(granges(GRanges(REs)), width=500, "center"), 
                   genome= "dm6", 
                   p.cutoff= 5e-4, 
                   bg= "even", 
                   out= "scores")
counts <- as.matrix(motifCounts(hit))
colnames(counts) <- name(vl_Dmel_motifs_DB$All_pwms_log_odds[sel])
counts <- as.data.table(counts)
names(counts) <- paste0("motif__", names(counts))

#-----------------------------------------#
# save
#-----------------------------------------#
REs <- cbind(REs, counts)
saveRDS(REs, 
        "Rdata/ED_REs_features_final.rds")






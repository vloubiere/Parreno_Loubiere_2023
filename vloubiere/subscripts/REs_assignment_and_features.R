setwd("D:/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

peaks <- fread("db/peaks/ATAC_seq_peaks_ED.txt")

genes <- import("D:/_R_data/genomes/dm6/dmel-all-r6.36.gtf")
genes <- as.data.table(genes)

#-----------------------------------------#
# Assemble RE table containing TSSs and non-TSS ATAC-Seq peaks
#-----------------------------------------#
TSS <- rbind(genes[type=="mRNA", .(seqnames= paste0("chr", seqnames), 
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
REs <- na.omit(REs)

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
                                 "D:/_R_data/projects/public_data/dm6/bw/merged_tracks/ATAC_ED_merged.bw")$score
REs$PH_score <- vl_bw_coverage(ext,
                               "D:/_R_data/projects/public_data/dm6/bw/SA_2020/PH_ED_merge.bw")$score
REs$SUZ12_score <- vl_bw_coverage(ext,
                                  "D:/_R_data/projects/public_data/dm6/bw/SA_2020/SUZ12_ED_merge.bw")$score
REs$H3K4me1_score <- vl_bw_coverage(ext,
                                    "D:/_R_data/projects/public_data/dm6/bw/SA_2020/H3K4me1_ED_merge.bw")$score
REs$H3K4me3_score <- vl_bw_coverage(ext,
                                    "D:/_R_data/projects/public_data/dm6/bw/SA_2020/H3K4me3_ED_merge.bw")$score
REs$H3K27ac_score <- vl_bw_coverage(ext,
                                    "D:/_R_data/projects/public_data/dm6/bw/SA_2020/H3K27Ac_ED_merge.bw")$score
REs$H3K27me3_score <- vl_bw_coverage(ext,
                                     "D:/_R_data/projects/public_data/dm6/bw/SA_2020/H3K27me3_ED_merge.bw")$score
REs$H3K27me2_score <- vl_bw_coverage(ext,
                                     "D:/_R_data/projects/public_data/dm6/bw/SA_2020/H3K27me2_ED_merge.bw")$score
REs$H2AK118Ub_score <- vl_bw_coverage(ext,
                                      "D:/_R_data/projects/public_data/dm6/bw/SA_2020/H2AK118Ub_ED_merge.bw")$score

#-----------------------------------------#
# CLEAN
#-----------------------------------------#
cols <- grep("score", colnames(REs), value = T)
REs <- REs[REs[, !is.na(rowSums(.SD)), .SDcols= cols]]
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
boxplot(mat)

pdf("pdf/clustering_REs.pdf")
cl <- vl_heatmap(mat, 
                 kmeans_k = 4, 
                 col = c("royalblue2", "yellow"))
dev.off()

cl <- cl[, switch(rcl, 
                  '1'="Enhancer", 
                  '2'="Polycomb", 
                  '3'="Null", 
                  '4'="Tss"), .(RE_ID= as.numeric(row), rcl)][order(RE_ID)]
REs[cl, HTM_cl:= i.V1, on= "RE_ID"]
setcolorder(REs, c("seqnames", "start", "end", "strand", "RE_ID", "HTM_cl"))

saveRDS(REs, 
        "Rdata/ED_REs_features_final.rds")






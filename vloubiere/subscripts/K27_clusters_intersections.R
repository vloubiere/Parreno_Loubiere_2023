setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(vlfunctions)
require(GenomicRanges)

# Import genes
gtf <- import("../../genomes/dm6/dmel-all-r6.36.gtf")
seqlevelsStyle(gtf) <- "UCSC"
gtf <- as.data.table(gtf)
genes_coor <- unique(gtf[type=="gene", .(seqnames, start, end, strand, FBgn= gene_id, gene_symbol)])

# Resize to extended promoters
genes_coor[strand=="+", c("start", "end"):= .(start-5000, start+2000)]
genes_coor[strand=="-", c("start", "end"):= .(end-2000, end+5000)]

# Import and format FC table
FC <- as.data.table(readRDS("Rdata/final_FC_table.rds"))
FC <- FC[grepl("PH.*TS", cdition)]
# FC[padj<0.05 & abs(log2FoldChange)>1, change:= ifelse(log2FoldChange>0, "up", "down")]
FC[padj<0.001, change:= ifelse(log2FoldChange>0, "up", "down")]
FC[is.na(change), change:= "unaffected"]
FC <- dcast(FC, FBgn~cdition, value.var= "change")
# Add transcriptome FC to genes
cols <- names(FC)[-1]
genes_coor[, (cols):= FC[.BY, .SD, .SDcols= cols, on= "FBgn"], FBgn]
genes_coor <- na.omit(genes_coor)

# Add K27 changes to genes
K27 <- fread("Rdata/K27_cutnrun_FC.txt")

cols <- grep("H3K27me3", names(K27), value = T)
K27me3 <- K27[region_type=="K27me3"]
K27me3[, (cols):= lapply(.SD, function(x) ifelse(abs(x)>1, ifelse(x>1, "up", "down"), "unaffected")), .SDcols= cols]
genes_coor <- cbind(genes_coor,
                    K27me3[genes_coor, {
                      # Identify max overlap for extraction
                      if(.N>0)
                      {
                        idx <- which.max(min(c(end, i.end))-max(c(start, i.start)))
                        .(H3K27me3= "bound",
                          H3K27me3_PH29_log2FC= H3K27me3_PH29_log2FC[idx], 
                          H3K27me3_PHD11_log2FC= H3K27me3_PHD11_log2FC[idx], 
                          H3K27me3_PHD9_log2FC= H3K27me3_PHD9_log2FC[idx])
                      }else
                      {
                        .(H3K27me3= "unbound",
                          H3K27me3_PH29_log2FC= "unaffected", 
                          H3K27me3_PHD11_log2FC= "unaffected", 
                          H3K27me3_PHD9_log2FC= "unaffected")
                      }
                    }, .EACHI, on= c("seqnames", "start<end", "end>start")][, .(H3K27me3, H3K27me3_PH29_log2FC, H3K27me3_PHD11_log2FC, H3K27me3_PHD9_log2FC)])

cols <- grep("H3K27Ac", names(K27), value = T)
K27Ac <- K27[region_type=="K27Ac"]
K27Ac[, (cols):= lapply(.SD, function(x) ifelse(abs(x)>1, ifelse(x>1, "up", "down"), "unaffected")), .SDcols= cols]
genes_coor <- cbind(genes_coor,
                    K27Ac[genes_coor, {
                      if(.N>0)
                      {
                        idx <- which.max(min(c(end, i.end))-max(c(start, i.start)))
                        .(H3K27Ac= "bound",
                          H3K27Ac_PH29_log2FC= H3K27Ac_PH29_log2FC[idx], 
                          H3K27Ac_PHD11_log2FC= H3K27Ac_PHD11_log2FC[idx], 
                          H3K27Ac_PHD9_log2FC= H3K27Ac_PHD9_log2FC[idx])
                      }else
                      {
                        .(H3K27Ac= "unbound",
                          H3K27Ac_PH29_log2FC= "unaffected", 
                          H3K27Ac_PHD11_log2FC= "unaffected", 
                          H3K27Ac_PHD9_log2FC= "unaffected")
                      }
                    }, .EACHI, on= c("seqnames", "start<end", "end>start")][, .(H3K27Ac, H3K27Ac_PH29_log2FC, H3K27Ac_PHD11_log2FC, H3K27Ac_PHD9_log2FC)])

# Add PRC1 binding to genes
PRC1 <- get(load("external_data/SA2020_cl.list"))
genes_coor[, PRC1:= ifelse(gene_symbol %in% unlist(PRC1$genes), "bound", "unbound")]

#-----------------------------#
# Test for associations
#-----------------------------#
comp <- data.table(cdition_1= c("PH29_TS", "PHD9_TS", "PHD11_TS"),
                   cdition_2= list(c("H3K27me3_PH29_log2FC", "H3K27Ac_PH29_log2FC", "PRC1"),
                                   c("H3K27me3_PHD9_log2FC", "H3K27Ac_PHD9_log2FC", "PRC1"),
                                   c("H3K27me3_PHD11_log2FC", "H3K27Ac_PHD11_log2FC", "PRC1")))
comp <- comp[, .(cdition_2= unlist(cdition_2)), cdition_1]
comp <- comp[, CJ(var1= genes_coor$PH29_TS,
                  var2= if(grepl("log2FC", cdition_2)) genes_coor$PH29_TS else c("bound", "unbound"), 
                  unique = T), .(cdition_1, cdition_2)]
comp <- comp[, .(subset= c("all_genes", "K27me3_genes", "K27Ac_genes")), (comp)]

# Compute fisher tests
comp[, c("estimate", "pval"):= {
  if(subset=="all_genes")
    .genes_coor <- genes_coor else if(subset=="K27me3_genes")
      .genes_coor <- genes_coor[H3K27me3=="bound"] else if(subset=="K27Ac_genes")
        .genes_coor <- genes_coor[H3K27Ac=="bound"]
  tab <- table(.genes_coor[[cdition_1]]==var1,
               .genes_coor[[cdition_2]]==var2)
  if(identical(dim(tab), c(2L, 2L)))
  {
    .c <- fisher.test(tab)
    .(.c$estimate , .c$p.value)
  }else
    .(as.numeric(NA), as.numeric(NA))
}, (comp)]
comp[, "-log10_padj":= -log10(p.adjust(pval, method= "fdr"))]
comp[, "log2_OR":= log2(estimate)]

# PLOT
pdf("pdf/cutnrun/fisher_associations_K27_vars.pdf", width = 13, height = 11)
par(mfrow=c(3, 3))
cditions <- data.table(cdition= c("PH29", "PHD11", "PHD9"))
cditions <- cditions[, .(c_subset= c("all_genes", "K27me3_genes", "K27Ac_genes")), (cditions)]
cditions[, {
  .c <- comp[grepl(cdition, cdition_1) & subset==c_subset]
  pl1 <- dcast(.c, 
               cdition_1+var1~cdition_2+var2, 
               value.var = "log2_OR")
  pl2 <- dcast(.c, 
               cdition_1+var1~cdition_2+var2, 
               value.var = "-log10_padj")
  pl1 <- t(as.matrix(pl1[, -1], 1))
  pl2 <- t(as.matrix(pl2[, -1], 1))
  pl1[pl1<0 | pl2<3] <- NA
  rownames(pl1) <- gsub(paste0("_", cdition, "_log2FC_"), " ", rownames(pl1))
  vl_balloons_plot(pl1,
                   pl2, 
                   x_breaks = c(0, 3), 
                   col= c("blue", "red"), 
                   main= paste(cdition, c_subset), 
                   balloon_size_legend = "log2(Odd ratio)",
                   balloon_col_legend = "-log10(padj)")
}, cditions]
dev.off()

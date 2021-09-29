gtf <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
ph29_DNA <- fread("db/DNA_analysis_novogene/customized_analysis/fisherresult/ph29_VS_ph29_t5_qvalue.txt")
ph29_RNA <- fread("db/FC_tables/RNA_epiCancer_PH29_T_FACS_vs_PH29_FC.txt")

gene_body <- GRanges(gtf[type=="gene" & gene_id %in% ph29_RNA[!is.na(pvalue), V1], .(seqnames, start, end, strand, gene_id, gene_symbol)])
ext_gene_body <- resize(gene_body, width(gene_body)+2000, "end")
promoter <- resize(gene_body, 1, fix = "start")
ext_promoter <- resize(promoter, 2000, "end")

res <- as.data.table(gene_body)
res[ph29_RNA, RNA_basemean:= log2(baseMean), on="gene_id==V1"]
res[ph29_RNA, RNA_log2FC:= log2FoldChange, on="gene_id==V1"]
res[ph29_RNA, RNA_log10pval:= -log10(i.pvalue), on="gene_id==V1"]
res$DNA_ext_gene_body_log10pval <- ph29_DNA[as.data.table(ext_gene_body), -log10(mean(fisher)), .EACHI, on= c("CHROM==seqnames", "POS>=start", "POS<=end")]$V1
res$DNA_ext_promoter_log10pval <- ph29_DNA[as.data.table(ext_promoter), -log10(mean(fisher)), .EACHI, on= c("CHROM==seqnames", "POS>=start", "POS<=end")]$V1

pdf("pdf/association_mutations_FC.pdf", height = 3.6, width = 10)
par(mfrow= c(1,3), las= 1)
# extended gene body
plot(res$RNA_log10pval, 
     res$DNA_ext_gene_body_log10pval, 
     xlim= c(0,100),
     xlab= "RNA -log10(pval)",
     ylab= "DNA -log10(pval) extended gene body",
     main= "ph29 vs ph29_t5")
pval <- fisher.test(table(res$RNA_log10pval>10, 
                          res$DNA_ext_gene_body_log10pval>10))$p.value
legend("topright", 
       legend = paste0("fisher=", formatC(pval, format = "e", digits = 2)),
       bty= "n")
# extended promoter
plot(res$RNA_log10pval, 
     res$DNA_ext_promoter_log10pval, 
     xlim= c(0,100),
     xlab= "RNA -log10(pval)",
     ylab= "DNA -log10(pval) extended promoter",
     main= "ph29 vs ph29_t5")
pval <- fisher.test(table(res$RNA_log10pval>10, 
                          res$DNA_ext_promoter_log10pval>10))$p.value
legend("topright", 
       legend = paste0("fisher=", formatC(pval, format = "e", digits = 2)),
       bty= "n")
# log2FC
plot(res$RNA_log2FC, 
     res$DNA_ext_gene_body_log10pval, 
     xlab= "RNA log2FoldChange",
     ylab= "DNA -log10(pval) extended gene body",
     main= "ph29 vs ph29_t5")
pval <- fisher.test(table(res$RNA_log2FC>2, 
                          res$DNA_ext_gene_body_log10pval>10))$p.value
legend("topright", 
       legend = paste0("fisher=", formatC(pval, format = "e", digits = 2)),
       bty= "n")
pval <- fisher.test(table(res$RNA_log2FC<(-2), 
                          res$DNA_ext_gene_body_log10pval>10))$p.value
legend("topleft", 
       legend = paste0("fisher=", formatC(pval, format = "e", digits = 2)),
       bty= "n")
dev.off()

man <- fread("db/DNA_analysis_novogene/customized_analysis/qvalueplot/ph29_VS_ph29_t5_meanq.txt")
man <- man[order(CHR, START)]
genome <- as.data.table(GRanges(GenomeInfoDb::seqinfo(BSgenome.Dmelanogaster.UCSC.dm6)))
genome[, seqnames:= gsub("chr", "", seqnames)]
genome <- genome[seqnames %in% man$CHR]
man[genome, x:= round((rowMeans(.SD)/i.width)*i.width), on= "CHR==seqnames", .SDcols= c("START", "END")]

ph29_RNA[as.data.table(promoter), c("seqnames", "start"):= .(i.seqnames, i.start), on= "V1==gene_id"]
ph29_RNA <- ph29_RNA[order(seqnames, start)]
ph29_RNA[genome, x:= round((start/i.width)*i.width), on= "seqnames"]

pdf("pdf/manhattan_plot_mutations.pdf", width = 10, height = 8)
par(mfrow= c(2,1),
    las= 1,
    mar= c(4,4,2,1))
for(chr in unique(man$CHR))
{
  plot(man[CHR==chr, x], 
       man[CHR==chr, `mean -log10(qvalue)/10Kb`], 
       type= "l",
       ylab= "DNA -log10(mean qvalue ph29T5 vs ph29)",
       main= chr,
       xlab= "Genomic coordinates")
  plot(ph29_RNA[seqnames==chr, x], 
       ph29_RNA[seqnames==chr, ifelse(log2FoldChange>0, -log10(padj), log10(padj))], 
       type= "l",
       ylab= "RNA -log10(mean padj ph29T5 vs ph29)",
       xlab= "Genomic coordinates")
}
dev.off()



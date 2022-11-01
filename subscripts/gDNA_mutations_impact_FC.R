setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(data.table)

FC <- fread("Rdata/final_gene_features_table.txt")
FC <- melt(FC, id.vars = "FBgn", patterns("diff"))
FC[, cdition:= gsub("^diff_", "", variable)]
FC <- FC[cdition!="PH18"]

mut <- fread("Rdata/mutated_gene_counts.txt")
mut <- melt(mut, 
            id.vars = c("FBgn", "symbol"), 
            measure.vars = patterns("count"= "^InDel|^SNP"))
mut[, cdition:= tstrsplit(variable, "_", keep= 2)]
mut <- mut[, .(value= sum(value)), .(FBgn, symbol, cdition)]
mut <- mut[value>0]

FC[, mutated:= FALSE]
FC[mut, mutated:= TRUE, on= c("FBgn", "cdition")]

FC[, c("OR_up", "pval_up"):= 
     fisher.test(value=="up", 
                 mutated, 
                 alternative = "greater")[c("estimate", "p.value")], cdition]
FC[, c("OR_down", "pval_down"):= 
     fisher.test(value=="down", 
                 mutated, 
                 alternative = "greater")[c("estimate", "p.value")], cdition]

res <- unique(FC[, .(variable, OR_up, pval_up, OR_down, pval_down)])
cols <- grep("pval", names(res), value = T)
res[, (cols):= lapply(.SD, function(x) formatC(x, format= "e", digits = 1)), .SDcols= cols]
cols <- grep("OR", names(res), value = T)
res[, (cols):= lapply(.SD, function(x) round(x, 2)), .SDcols= cols]

pdf("pdf/gDNA_mutations_impact_FC.pdf", width = 10, height = 5)
vl_plot_table(FC[, .N, keyby= .(variable, value, mutated)])
vl_plot_table(res)
dev.off()

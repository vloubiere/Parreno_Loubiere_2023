
GO <- get_ontology("D:/_R_data/public_data/dm6/GO/go-basic.obo")
clean <- data.table(id= GO$id, name= GO$name)
clean <- clean[grepl("^GO", id)]

current <- fread("D:/_R_data/public_data/dm6/GO/gene_association_FB2021_01.fb", skip = 5)
current[clean, name:= i.name, on= "V5==id"]
current[, total:= .N, V5]
sub <- unique(current[V3 %in% c("l(2)gl", "dlg1", "crb", "scrib", "Moe"), .(V3, V5, name, total)])
sub <- sub[, .(genes= .(unique(V3))), .(V5, name, total)]
sub[, N:= lengths(genes)]
sub[order(-N)][1:20]

res <- unique(current[name== "apicolateral plasma membrane", .(symbol= V3, FBgn= V2)])
fwrite(res, "JP_data_RNA_PolyComb_01.02/polarity_genes.txt", sep= "\t")

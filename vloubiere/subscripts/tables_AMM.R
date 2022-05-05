setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")

#----------------------------------------------------------#
# Make tables AMM
#----------------------------------------------------------#
tables <- list.files("db/FC_tables/RNA/", "epiCancer", full.names = T)
tables <- rbindlist(lapply(setNames(tables, gsub("epiCancer_ED_|_RNA_gDNA_RNA|.txt|RNA_","",basename(tables))), fread), 
                    idcol = T)
tables <- dcast(tables[, .(.id, FBgn, log2FoldChange, padj)], 
                FBgn~.id, 
                value.var = list("log2FoldChange", "padj"))
tables[as.data.table(rtracklayer::import("../../genomes/dm6/dmel-all-r6.36.gtf")), symbol:= gene_symbol, on= "FBgn==gene_id", mult="first"]
setcolorder(tables, 
            c("FBgn", 
              "symbol", 
              c(matrix(c(grep("log2FoldChange", names(tables), value= T),
                         grep("padj", names(tables), value= T)), 
                       nrow= 2, 
                       byrow= T))))

# Add clusters
cl <- readRDS("Rdata/clustering_RNA.rds")
tables[cl$allograft$rows, allograft_cl:= i.cl, on="FBgn"]
tables[cl$cutnrun$rows, cutnrun_genotype_cl:= i.cl, on="FBgn"]

# Add PRC1 binding
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}
PRC1 <- loadRData("external_data/SA2020_cl.list")
PRC1 <- rbindlist(lapply(PRC1$genes, function(x) data.table(symbol= x)), idcol = "PRC1_cluster")
tables[PRC1, PRC1_cluster:= i.PRC1_cluster, on= "symbol"]

# Order and save
setcolorder(tables, 
            c("FBgn", "symbol", "allograft_cl", "cutnrun_genotype_cl", "PRC1_cluster"))
fwrite(tables, 
       "db/FC_tables/RNA_table_AMM.txt",
       col.names = T, 
       row.names = F, 
       sep= "\t",
       quote= F, 
       na = NA)

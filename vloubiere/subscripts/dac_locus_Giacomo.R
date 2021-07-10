files <- as.list(c(list.files("db/FC_tables/", full.names = T), 
                   list.files("external_data/", "Bilder", full.names = T)))
dat <- lapply(files, function(x) 
{
  if(grepl(".txt$", x))
    x <- fread(x)
  else if(grepl(".xlsx", x))
  {
    x <- as.data.table(read_xlsx(x))
    colnames(x)[7] <- "log2FoldChange"
  }
  colnames(x)[1] <- "FBgn"
  return(x)
})
names(dat) <- files
dat <- rbindlist(dat, 
                 idcol = "cdition", 
                 fill= T)
dat[, cdition:= gsub("_ED|ED|_FC.txt$", "", basename(cdition)), cdition]
dat[, cdition:= gsub("2018_|RNA_", "", cdition), cdition]
dat[, cdition:= gsub("phRNAi_SA2020", "25C", cdition), cdition]
dat[, cdition:= gsub("mutants_SA2020", "CpGmut", cdition), cdition]
dat[, cdition:= gsub("^(.*)_transcriptome_Bilder.*", "Bilder_\\1", cdition), cdition]
dat[, log2FoldChange:= as.numeric(log2FoldChange)]

gtf <- as.data.table(import("../../genomes/dm6/dmel-all-r6.36.gtf"))
dat[gtf[type=="gene"], symbol:= i.gene_symbol, on= "FBgn==gene_id"]

res <- dat[symbol %in% c("CG5888", "Idgf1", "dac"), .(symbol, cdition, log2FoldChange)]
res <- as.matrix(dcast(res, symbol~cdition, value.var = "log2FoldChange"), 1)

pdf("pdf/dac_locus_giacomo.pdf", width = 10, height = 6)
par(mar= c(20,5,3,3))
vl_heatmap(res, 
           display_numbers = T, 
           cluster_cols = F)
dev.off()

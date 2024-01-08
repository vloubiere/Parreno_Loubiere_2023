setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# FC
dat <- data.table(file= list.files("db/FC_tables/RNA/", "rescue", full.names = T))
dat[, name:= gsub("epiCancer_rescue_|.txt", "", basename(file))]
dat <- dat[, fread(file), name]
gtf <- rtracklayer::import("../../genomes/Drosophila_melanogaster/flybase/dm6/dmel-all-r6.36.gtf")
gtf <- as.data.table(gtf)
dat[gtf, symbol:= i.gene_symbol, on= "FBgn==gene_id"]
dat[, diff:= fcase(log2FoldChange>1 & padj<0.05, "up",
                   log2FoldChange>1 & padj<0.05, "down",
                   default = "unaffected")]
setcolorder(dat, c("name", "FBgn", "symbol"))
dat <- split(dat[, !"name"], dat$name)

#Add FPKM
FPKM <- readRDS("db/dds/RNA/epiCancer_rescue_dds.rds")
FPKM <- as.data.table(DESeq2::fpkm(FPKM), keep.rownames = "FBgn")
FPKM <- melt(FPKM, id.vars = "FBgn")
FPKM[, c("variable", "rep"):= tstrsplit(variable, "__")]
FPKM <- FPKM[, .(value= mean(value)), .(FBgn, variable)]
FPKM <- dcast(FPKM, FBgn~variable, value.var = "value")
dat <- c(dat, list(FPKM= FPKM))

writexl::write_xlsx(dat,
                    "db/extended_data_table/Extended_Data_Table_4.xlsx")

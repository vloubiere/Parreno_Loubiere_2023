setwd("D:/_R_data/projects/epigenetic_cancer/")
require(data.table)
require(vlfunctions)

dat <- data.table(file= list.files("db/counts", 
                                   recursive = T, 
                                   full.names = T))
meta <- fread("Rdata/raw_metadata_final.txt")
meta[, counts:= paste0("db/counts/", project, "/", prefix_output, "_counts.rds")]
dat[meta, c("DESeq2_object", "cdition"):= .(i.DESeq2_object, 
                                            i.prefix_output), on= "file==counts"]
dat <- dat[, .(DESeq2_object= unlist(tstrsplit(DESeq2_object, ";"))), setdiff(colnames(dat), "DESeq2_object")]
dat <- dat[, {
        .c <- readRDS(file)$counts
        data.table(rn= rownames(.c), 
                   counts= .c[,1])
}, (dat)]

sdat <- split(dat, dat$DESeq2_object)
sdat <- lapply(sdat, function(x)
{
        .c <- dcast(x, rn~cdition, 
                    value.var = "counts")
        return(.c[rowSums(.c[, -1])>0])
})
        
dir.create("pdf/PCC_heatmaps", 
           showWarnings = F)

for(i in seq(sdat))
{
        .c <- cor(as.matrix(sdat[[i]], 1))
        pdf(paste0("pdf/PCC_heatmaps/", names(sdat)[i], ".pdf"),
            width = 10,
            height = 10)
        vl_heatmap(as.matrix(.c), 
                   display_numbers = T)
        dev.off()
}

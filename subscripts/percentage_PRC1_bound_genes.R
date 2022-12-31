# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import data
dat <- fread("Rdata/final_gene_features_table.txt")
dat[is.na(cl), cl:= 0]
res <- list()
for(ccl in 0:6)
{
  .c <- rbind(dat[cl==0, .(PRC1_bound)],
              dat[cl==ccl, .(PRC1_bound)])
  res[[ccl+1]] <- fisher.test(c(rep(F, sum(dat$cl==0)), rep(T, sum(dat$cl==ccl))),
                              .c$PRC1_bound)[c("estimate", "p.value")]
}
names(res) <- paste0(0:6)
pl <- rbindlist(res, idcol = T) 
dat[, cl:= as.character(cl)]
dat <- dat[pl, on= "cl==.id"]
dat[, total_cluster:= .N, cl]
pl <- dat[(PRC1_bound), .N, .(cl, estimate, p.value, total_cluster, K27me3_bound)]
pl[, perc:= N/total_cluster*100, .(cl, total_cluster)]


pdf("pdf/Extended_data_5_cluster_PRC1_K27me3_binding.pdf",
    width = 2.5,
    height = 2.5)
par(mar= c(3,3,1,0.5),
    mgp= c(2,0.5,0),
    las= 1,
    tcl= -0.2)
bar <- barplot(perc~K27me3_bound+cl, 
        pl,
        border= NA,
        ylab= "PRC1 bound genes (%)",
        ylim= c(0, 40), 
        beside= F, 
        col= c("lightgrey", "cornflowerblue"),
        xlab= NA,
        xaxt= "n")
vl_tilt_xaxis(x= bar, 
              labels= c("Unaffected", paste0("Cluster", 1:6)),
              srt= 30)
pl[, vl_plot_pval_text(bar[.GRP], 
                       sum(perc), 
                       p.value, 
                       stars_only = T), cl]
legend("topleft",
       fill= c("lightgrey",
               "cornflowerblue"),
       legend= c("H3K27me3-",
                 "H3K27me3+"), 
       bty= "n",
       cex= 0.6)
dev.off()
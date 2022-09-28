setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
require(vlfunctions)
require(data.table)

# Import recovery classes and transcriptome clusters
dat <- fread("Rdata/final_gene_features_table.txt")
dat[, chromhmm:= factor(chromhmm)]

###########################################
# PLOT
###########################################
pdf("pdf/recovery_chromhmm_class.pdf", 
    width = 5,
    height = 5)
par(las= 2,
    mar= c(8,10,8,10),
    tcl= -0.2,
    mgp= c(2,0.5,0))
counts <- do.call(cbind,
                     list(all_genes= table(dat$chromhmm),
                          recovery= table(dat[recovery=="Recovery", chromhmm]),
                          NoRecovery= table(dat[recovery=="noRecovery", chromhmm])))
perc <- apply(counts, 2, function(x) x/sum(x)*100)
Cc <- c("tomato", "brown", "forestgreen", "grey", "cornflowerblue")
Cc <- adjustcolor(Cc, 0.5)
bar <- barplot(perc, 
               col= Cc,
               ylab= "Percentage of genes",
               main= "chromHMM SA 2020")
counts[perc<5] <- NA
text(rep(bar, each= nrow(perc)),
     apply(perc, 2, cumsum)-perc/2,
     counts,
     cex= 0.5)
legend(par("usr")[2],
       par("usr")[4],
       legend= rownames(perc),
       fill= Cc,
       xpd= T, 
       bty= "n")
dev.off()
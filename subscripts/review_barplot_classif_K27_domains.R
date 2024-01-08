setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import data ----
dat <- readRDS("Rdata/K27_domains_classif_reversible_irreversible.rds")
dat <- dat[irreversible+reversible>0]
dat[rev_symbol!="" & irrev_symbol!="", irrev_symbol:= paste0(irrev_symbol, ", ")]
dat[, rev_symbol:= gsub("lncRNA:|asRNA:", "", rev_symbol)]
dat[, irrev_symbol:= gsub("lncRNA:|asRNA:", "", irrev_symbol)]

# Plot ----
Cc <- c("rosybrown1", "palegreen3", "lightgrey")

pdf("pdf/review_K27me3_domains_reversible_irreversible.pdf", 3, 4.5)
par(mai= c(.5,2,0.2,.2),
    tcl= -0.1,
    mgp= c(.5,0.35,0),
    las= 1,
    cex.lab= 8/12,
    cex.axis= 7/12)
bar <- barplot(t(dat[, irreversible:other]),
               border= NA,
               col= Cc,
               horiz= T,
               xlab= "Number of genes",
               xaxt= "n")
text(par("usr")[1],
     bar,
     dat$rev_symbol,
     offset= 0,
     pos= 2,
     xpd= T,
     cex= 5/12,
     col= Cc[2])
text(par("usr")[1]-strwidth(dat$rev_symbol, cex= 5/12),
     bar,
     dat$irrev_symbol,
     offset= 0, 
     pos= 2,
     xpd= T,
     cex= 5/12,
     col= Cc[1])
axis(1, padj= -1.5)
legend(grconvertX(0.1, "nfc", "user"),
       par("usr")[4],
       fill= Cc,
       legend= c("Irreversible", "Reversible", "Other"),
       cex= 8/12,
       border= NA,
       bty= "n",
       xpd= NA)
dev.off()
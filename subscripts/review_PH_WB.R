setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")

# Import
dat <- fread("db/WB/ph_WB_fig1_new.txt")
# Normalize
dat[, norm:= `PH/TUB`/mean(`PH/TUB`[cdition=="No ph-KD"])]
dat <- dat[, .(mean= mean(norm),
               sd= sd(norm),
               value= .(norm)), cdition]


pdf("pdf/review_WB_quantif_barplot.pdf", 3, 3)
vl_par(xpd= T)
bar <- barplot(dat$mean,
               ylab= "Normalized PH/TUB ratio")
points(rep(bar, each= 3),
       unlist(dat$value),
       pch= 16,
       col= adjustcolor("black", .6))
arrows(bar,
       dat$mean,
       bar,
       dat[,mean-sd],
       angle = 90,
       length = .05)
arrows(bar,
       dat$mean,
       bar,
       dat[,mean+sd],
       angle = 90,
       length = .05)
segments(bar[1], 1.1, bar[2], 1.1)
vl_plot_pval_text(mean(bar[c(1, 2)]),
                  1.1,
                  t.test(dat$value[[1]], dat$value[[2]])$p.value,
                  stars.only = T)
segments(bar[1], 1.2, bar[3], 1.2)
vl_plot_pval_text(mean(bar[c(1, 3)]),
                  1.2,
                  t.test(dat$value[[1]], dat$value[[3]])$p.value,
                  stars.only = T)
segments(bar[1], 1.3, bar[4], 1.3)
vl_plot_pval_text(mean(bar[c(1, 4)]),
                  1.3,
                  t.test(dat$value[[1]], dat$value[[4]])$p.value,
                  stars.only = T)
dev.off()
require(data.table)

dat <- fread("db/WB/ph_WB_fig1.txt", fill= T)
pl <- dat[, .(mean= mean(`ph/tub_ratio`), sd= sd(`ph/tub_ratio`)), Cdition]

pdf("pdf/Figure_1_WB_PH.pdf", 4, 3)
par(mgp= c(2, 0.5, 0),
    tcl= -0.2,
    las= 1,
    mar= c(7,5,5,2),
    xpd= NA)
bar <- barplot(pl$mean,
               border= NA,
               ylab= "PH/TUB ratio",
               width= 0.5, space= 1)
vl_tilt_xaxis(bar, labels= gsub("_", " ", pl$Cdition))
arrows(bar, pl$mean, bar, pl[, mean+sd], angle = 90, length = 0.1)
arrows(bar, pl$mean, bar, pl[, mean-sd], angle = 90, length = 0.1)

# Test constant
test <- t.test(dat[Cdition=="No_ph-KD", `ph/tub_ratio`], 
               dat[Cdition=="Constant_ph-KD", `ph/tub_ratio`])$p.value
segments(bar[1], 4.2, bar[2], 4.2)
vl_plot_pval_text(mean(bar[1:2]), 4.2, test, pos= 2, stars_only = T)

# Test d9
test <- t.test(dat[Cdition=="No_ph-KD", `ph/tub_ratio`], 
               dat[Cdition=="Transient_ph-KD_d9", `ph/tub_ratio`])$p.value
segments(bar[1], 4.8, bar[3], 4.8)
vl_plot_pval_text(mean(bar[1:3]), 4.9, test, pos= 2, stars_only = T)

# Test d11
test <- t.test(dat[Cdition=="No_ph-KD", `ph/tub_ratio`], 
               dat[Cdition=="Transient_ph-KD_d11", `ph/tub_ratio`])$p.value
segments(bar[1], 5.4, bar[4], 5.4)
vl_plot_pval_text(mean(bar[1:4]), 5.5, test, pos= 2, stars_only = T)

dev.off()

setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(ggalluvial)

# Import metadata ----
dat <- readRDS("Rdata/final_ATAC_rescue_table.rds")

# Only peaks that are decently strong and show a difference in at least one of the conditions ----
dat <- dat[log10(baseMean)>1.25 & (diff_gfp_ph!="unaffected" | diff_zfh1_ph!="unaffected")]

# Prepare for plotting ----
pl <- na.omit(dat[, .(diff_gfp_ph, diff_zfh1_ph)])
pl <- pl[, .(freq= .N, subject= .GRP), (pl)]
pl <- melt(pl,
            measure.vars = patterns("diff"))
pl[, value:= factor(value, c("up", "unaffected", "down"))]

 # Plot ----
pdf("pdf/review_alluvial_plot_rescue_ATAC.pdf",
    width = 3.5,
    height = 3.5)
ggplot(pl,
       aes(x = variable, stratum = value, alluvium = subject,
           y = freq,
           fill = value, label = value)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = c(up= "tomato", down= "cornflowerblue", unaffected= "grey80"))+
  geom_flow(alpha = 0.5, width= 1/1.5) +
  geom_stratum(alpha = .7, colour= NA, width = 1/1.5) +
  geom_text(stat = "stratum", size = 3) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.2,1,0.2,1), "cm"),
        legend.position = "none")
dev.off()
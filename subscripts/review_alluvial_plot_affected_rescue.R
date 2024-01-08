setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(ggalluvial)

# Import data
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[, .(diff_GFP_PH, diff_STAT92E_PH, diff_ZFH1_PH)]
dat <- dat[apply(dat, 1, function(x) any(x!="unaffected"))]
dat <- na.omit(dat)
dat <- dat[, .(freq= .N, subject= .GRP), (dat)]
dat <- melt(dat,
            measure.vars = patterns("diff"))
dat[, value:= factor(value, c("up", "unaffected", "down"))]

# Plot ----
pdf("pdf/review_alluvial_plot_rescue_transcriptomes.pdf", 
    width = 3.5, 
    height = 3.5)
ggplot(dat,
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
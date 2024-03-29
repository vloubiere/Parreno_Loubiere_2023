setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(ggalluvial)

# Import metadata
dat <- readRDS("Rdata/final_gene_features_table.rds")
dat <- dat[, .(diff_PH18, diff_PH29, diff_PHD9, diff_PHD11)]
dat <- dat[apply(dat, 1, function(x) any(x!="unaffected"))]
setnames(dat,
         c("no ph-KD", "Constant ph-KD", "Transient ph-KD d9", "Transient ph-KD d11"))
dat <- na.omit(dat)
dat <- dat[, .(freq= .N, subject= .GRP), (dat)]
dat <- melt(dat,
            measure.vars = setdiff(names(dat), c("freq", "subject")))
dat[, value:= factor(value, c("up", "unaffected", "down"))]

# Plot ----
pdf("pdf/Figure_2_alluvial_plot_timecourse_transcriptome.pdf", 
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
  ylab("Number of genes\n(vs temperature-matched w-KD)")+
  theme(axis.line = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust= 1),
        axis.title.x = element_blank(),
        plot.margin = unit(c(0.2,1,0.2,1), "cm"),
        legend.position = "none")
dev.off()

# setwd("/mnt/d/_R_data/projects/epigenetic_cancer/")
setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(vlfunctions)
require(readxl)
require(ggalluvial)

# Import metadata
meta <- fread("Rdata/processed_metadata_RNA.txt")
dat <- meta[!is.na(FC_file) & system=="noGFP", fread(FC_file), .(FC_file, cdition)]
dat[, cdition:= switch(cdition, 
                       "PH18"= "no ph-KD",
                       "PH29"= "Constant ph-KD",
                       "PHD9"= "Transient ph-KD d9",
                       "PHD11"= "Transient ph-KD d11"), cdition]
dat[, cdition:= factor(cdition, 
                       levels= c("no ph-KD", "Constant ph-KD", "Transient ph-KD d9", "Transient ph-KD d11"))]
dat <- dcast(dat,
             FBgn~cdition,
             value.var = "diff")
dat <- dat[, .(freq= .N, subject= .GRP), `no ph-KD`:`Transient ph-KD d11`]
dat <- melt(dat, measure.vars = patterns("ph-KD"))
dat[, value:= factor(value, c("up", "unaffected", "down"))]

pdf("pdf/Figure_2_alluvial_plot_timecourse_transcriptome.pdf", 
    width = 5, 
    height = 5)
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

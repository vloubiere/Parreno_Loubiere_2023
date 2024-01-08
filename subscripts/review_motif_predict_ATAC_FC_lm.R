setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
require(ggrepel)

# Import data ----
dat <- readRDS("db/model/ATAC_FC_PH29_PHD11_motif_lm.rds")
dat[, Cc:= fcase(`t value.d11`>5 & `t value.29`<(-5), "tomato",
                 `t value.d11`>5, "red",
                 `t value.d11`<(-7) & `t value.29`<(-7), "cornflowerblue",
                 `t value.d11`<(-7), "lightblue",
                 `t value.29`>7, "khaki3",
                 `t value.29`<(-7), "cyan3",
                 default= "lightgrey")]
dat[, Cc:= factor(Cc)]
dat <- dat[`Pr(>|t|).29`<1e-5 | `Pr(>|t|).d11`<1e-5]

# ggrepl
p <- ggplot(dat, aes(x = `t value.d11`, y = `t value.29`, color= Cc)) +
  geom_point() +
  scale_color_manual(values = levels(dat$Cc))+
  geom_text_repel(aes(label = name), box.padding = 0.15, point.padding = 0.15, size= 1.75) +
  theme_minimal() +
  # xlab("index") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),  # Add axis lines
    axis.ticks = element_line(color = "black"),  # Add ticks
    axis.ticks.length = unit(0.2, "cm")         # Set the length of ticks
  )+
  theme(legend.position = "none")

# Plot ----
pdf("pdf/review_motif_predict_ATAC_FC.pdf", 3, 3)
p
dev.off()
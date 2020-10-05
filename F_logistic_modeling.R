setwd("/groups/stark/vloubiere/projects/epigenetic_cancer/")
sapply(list.files("/groups/stark/vloubiere/functions/", ".R", full.names = T), source)
options(datatable.print.topn= 1)
require(data.table)
require(kohonen)
require(caret)

som <- readRDS("Rdata/som_genes_FC_PH_vs_W.rds")
dat <- readRDS("Rdata/genes_features_final.rds")
setkeyv(dat, "FBgn")
dat <- dat[rownames(som$data[[1]])]
dat[, cl:= som$kcl[som$unit.classif]]

# cols <- grep("_enr_", colnames(dat), value = T)
cols <- grep("_count_", colnames(dat), value = T)
sub <- dat[, ..cols]
sub[, class:= ifelse(dat$cl==1, "yes", "no")]
mem <- data.table(colnames(sub[, !"class"]))
mem[, V2:= paste0("col_", .I)]
colnames(sub)[-ncol(sub)] <- mem$V2

trctrl <- trainControl(method = "repeatedcv", number = 2,repeats = 1, classProbs = T, savePredictions = T)
cmb <- CJ(colnames(sub[, !"class"]), colnames(sub[, !"class"]))
cmb <- cmb[V1!=V2]
cmb[, inter:= paste0(V1, "*", V2), (cmb)]

# .f <- as.formula(paste("class~", paste(colnames(sub[, !"class"]), collapse= "+")))
# .f <- as.formula(paste("class~", paste(cmb$inter, collapse= "+")))
.f <- as.formula(paste0("class~", paste(unique(cmb$V1), collapse= "+")))
.glm <- train(.f, data = sub, method = "glm", family="binomial",
              trControl= trctrl, preProcess = c("center", "scale"), tuneLength = 1)

test <- as.data.table(summary(.glm)$coefficients, keep.rownames = T)

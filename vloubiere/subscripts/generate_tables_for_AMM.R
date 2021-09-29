setwd("/mnt/d/_R_data/projects/epigenetic_cancer")
require(xlsx)
require(openxlsx) #recall the library

#--------------------#
# IMPORT FC
#--------------------#
df <- readRDS("Rdata/final_FC_table.rds")
df <- dcast(df, 
            FBgn+symbol~cdition, 
            value.var = list("log2FoldChange", "padj"))

#--------------------#
# IMPORT Clusters and add to DF
#--------------------#
dose_cl <- readRDS("Rdata/clustering_dose_transcriptomes.rds")
dose_cl <- unique(dose_cl[, .(rcl, FBgn)])
df[dose_cl, dose_cl:= i.rcl, on= "FBgn"]
epi_cl <- readRDS("Rdata/clustering_epiCancer_transcriptomes.rds")
epi_cl <- unique(epi_cl[, .(rcl, FBgn)])
df[epi_cl, epiCancer_cl:= i.rcl, on= "FBgn"]

names(df)[grep("^log2FoldChange", names(df))] <- gsub("^log2FoldChange_(.*)", "\\1_log2FC", names(df)[grep("^log2FoldChange", names(df))])
names(df)[grep("^padj", names(df))] <- gsub("^padj_(.*)", "\\1_padj", names(df)[grep("^padj", names(df))])

#--------------------#
# IMPORT gene features and add to data
#--------------------#
feat <- readRDS("Rdata/ED_REs_features_final.rds")
feat <- feat[, .(TSS_open_PRC1_total= paste0(length(which(grepl("TSS", type) & open)), "_", 
                                             length(which(grepl("TSS", type) & PRC1_bound)), "_", 
                                             length(grep("TSS", type))),
                 distalRE_open_PRC1_total= paste0(length(which(grepl("ENHANCER", type) & open)), "_" ,
                                                  length(which(grepl("ENHANCER", type) & PRC1_bound)), "_", 
                                                  length(grep("ENHANCER", type)))), FBgn]
df <- merge(df, 
            feat, 
            by= "FBgn")

#--------------------#
# WT FPKM
#--------------------#
fpkm <- fread("db/fpkms/RNA_epiCancer_ED_handDissect_Parreno.txt")
fpkm[, log2_fpkm_white:= apply(.SD, 1, mean, na.rm= T), .SDcols= patterns("^W")]
df[fpkm, log2_fpkm_white:= i.log2_fpkm_white, on= "FBgn"]

setcolorder(df, 
            c("FBgn",
              "symbol", 
              "dose_cl", 
              "epiCancer_cl", 
              "TSS_open_PRC1_total", 
              "distalRE_open_PRC1_total",
              "log2_fpkm_white"))

#--------------------#
# Format excel file
#--------------------#
# create a workbook
wb <- createWorkbook() 
#add a worksheet to the workbook
addWorksheet(wb, "Sheet", 
             gridLines = TRUE) 
# write my analysis into the worksheet of the workbook
writeData(wb, 
          "Sheet", 
          df, 
          keepNA = T) 
# Create a list of cell format depending on padj/FC
Cc <- colorRampPalette(c("cornflowerblue", "white", "tomato"))(7)
checks <- list(down3= list(test= "AND($col11<0,$col21<0.05)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[3])),
               down2= list(test= "AND($col11<0,$col21<0.01)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[2])),
               down1= list(test= "AND($col11<0,$col21<0.001)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[1])),
               down1= list(test= "AND($col11>0,$col21<0.05)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[5])),
               down2= list(test= "AND($col11>0,$col21<0.01)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[6])),
               down3= list(test= "AND($col11>0,$col21<0.001)",
                           style= createStyle(fontColour = "#000000", 
                                              bgFill = Cc[7])))
# Create canonical excel col names
columns <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
columns <- c(columns, paste0("A", columns))
columns <- c(columns, paste0("B", columns))
# account padj and FC to colour cells
for(icol in grep("_log2FC$", names(df)))
{
  for(j in seq(checks))
  {
    conditionalFormatting(wb, 
                          "Sheet", 
                          cols = icol,
                          rows = 1:nrow(df),
                          rule = gsub("col2", columns[icol+12], gsub("col1", columns[icol], checks[[j]]$test)), 
                          style = checks[[j]]$style,
                          type = "expression")
  }
}
dir.create("tables_AMM", 
           showWarnings = F)
saveWorkbook(wb, 
             "tables_AMM/final_table.xlsx", 
             overwrite = TRUE)

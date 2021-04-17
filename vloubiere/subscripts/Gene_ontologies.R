# Generate GO table
if(!file.exists("Rdata/fb_go_table.rds"))
{
  keep.genes <- fread("db/FC_tables/RNA_epiCancer_PH18_vs_Ez18_FC.txt")$V1
  go <- fread("../../genomes/dm6/FB2020_05_gene_association.fb", 
              skip= 5, 
              select = c(2,3,5))
  colnames(go) <- c("Fbgn", "symbol", "GO")
  go <- go[Fbgn %in% keep.genes]
  go_details <- get_ontology("../../genomes/dm6/FB2020_05_go-basic.obo", 
                             extract_tags = "everything")
  go_details <- data.table(GO= go_details$id, 
                           name= go_details$name, 
                           type= go_details$namespace, 
                           parents= go_details$parents, 
                           children= go_details$children)
  go <- go_details[go, , on= "GO"]
  saveRDS(go, "Rdata/fb_go_table.rds")
}

# Compute GO ewnrichments
go <- readRDS("Rdata/fb_go_table.rds")
dat <- readRDS("Rdata/final_clustering_transcriptomes.rds")
go[, genes_GO:= length(unique(Fbgn)), GO]
go[, genes_noGO:= length(unique(go$Fbgn))-genes_GO]
enrich <- dat[, {
  .c <- unique(FBgn)
  res <- go[.c, .(genes_GO, 
                  genes_noGO, 
                  cl_GO= length(unique(Fbgn))), 
            on= "Fbgn", 
            by= .(GO, name), 
            nomatch= NULL]
  res[, cl_noGO:= length(.c)-cl_GO]
  res <- unique(res)
  res[, fisher.test(data.table(c(cl_GO, cl_noGO), 
                               c(genes_GO, genes_noGO)))[c("estimate", "p.value")], (res)]
}, cl]
enrich[, p.adj:= p.adjust(p.value, method = "fdr")]
enrich[, log2FE:= log2(estimate)]

# Select signif ontologies
sel <- enrich[, any(p.adj<0.05 & log2FE>0), GO][(V1), GO]
dmat <- enrich[GO %in% sel]
setorderv(dmat, c("cl", "log2FE"), c(1, -1))
dmat[, x:= .GRP, .(cl, log2FE)]
dmat[, x:= min(.SD[p.adj<0.05, x]), name]
dmat[, x:= .GRP, keyby= x]
setorderv(dmat, "cl", order = -1)
dmat[, y:= .GRP, cl]

# PLOT
pdf("pdf/balloon_GO_clusters.pdf", width = 15)
par(mar= c(15,5,2,2), xaxs= "i")
Cc <- cut(dmat$p.adj, 
          c(-Inf, 0.00001, 0.001, 0.01, 0.05, Inf), 
          c("red", "tomato", "cornflowerblue", "blue", "black"))
plot.new()
abline(v= unique(dmat$x/max(dmat$x)), lwd= 0.5)
abline(h= unique(dmat$y/max(dmat$y)), lwd= 0.5)
points(x= dmat$x/max(dmat$x), 
       y= dmat$y/max(dmat$y), 
       cex= dmat$log2FE/2.5, 
       col= as.character(Cc), 
       pch= 16, 
       xpd= T)
axis(2, 
     at = unique(dmat[, .(cl, y/max(y))])$V2, 
     labels = unique(dmat[, .(cl, y)])$cl, 
     las= 1, 
     tick= 0, 
     lwd= 0)
par(cex= 0.7)
axis(1, 
     at = unique(dmat[, .(name, x/max(x))])$V2, 
     labels = unique(dmat[, .(name, x)])$name, 
     las= 2, 
     tick= 0, 
     lwd= 0)
dev.off()



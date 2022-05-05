require(qpdf)

qpdf::pdf_combine(list.files("pdf/Figures/", "_RNA.pdf$", full.names = T),
                  "pdf/Figures/all_figures_RNA.pdf")
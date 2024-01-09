# git_epiCancer

This repository contains all custom R (version 4.2.0) scripts supporting the conclusions of the study entitled "Transient loss of Polycomb components induces an epigenetic cancer fate".
Authors: Parreno V., Loubiere V., Schuettengruber B., Fritsch L., Rawal C. C., Erokhin M., Győrffy B., Normanno D., Di Stefano M., Moreaux J., Butova N. L, Chiolo I, Chetverina D., Martinez A-M & Cavalli G.
BioXriv doi: https://doi.org/10.1101/2023.01.04.522799.

The "main.R" file lists all the different scripts (see the "subscripts/" folder) that were generated, starting from the processing of raw sequencing reads and up to figure panels.
Script names should be self-explanatory, with additional comments linking them to the corresponding figure panel when relevant.

All custom functions that were generated were wrapped into a R package that was made publicly available: https://github.com/vloubiere/vlfunction/tree/nature_v2_revised.
The package can be installed using: install_github("https://github.com/vloubiere/vlfunction/tree/nature_v2_revised")

For any reasonable further request, please contact anne-marie.martinez@igh.cnrs.fr, giacomo.cavalli@igh.cnrs.fr, vincent.loubiere@imp.ac.at

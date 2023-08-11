# Code manuscript: Dual immune checkpoint blockade induces analogous alterations in the dysfunctional CD8+ T cell and activated Treg compartment

Code belonging to the manuscript "Dual immune checkpoint blockade induces analogous alterations in the dysfunctional CD8+ T cell and activated Treg compartment" (van der Leun et al. 2023). 

### Stucture code
1. Code related to the annotation of the immune cells using the **Metacell** R package (Baran et al. 2019) can be found in the **code/preprocessing** folder. Specifcally, `1.separate_objects.R`, `2.hn_full_dataset_merging.R` and `3.hn_metacell_generation_coley.R`.
2. Code related to reproducing the analysis and figures can be found in the **code/figures** folder, structured by individual figure.
3. Code related to the preprocessing of the Kurten et al. (Kürten et al. 2021) data can be found in the **code/preprocessing** folder, `1.kurten_metacell_generation.R`

### Data
- The processed single cell RNA sequencing data and clonotype information have been deposited in the NCBI GEO database (**GSE232240**). Count and annotation data are equal to the combined seurat objects in the preprocessing step in 1.separate_objects.R.
- The raw single cell RNA sequencing and TCR sequencing data have been deposited in the EGA database (**EGAS00001007367** & **EGAS00001007368**).

### Notes
- All code was written in R in Rstudio (v4.0.2). 

### References
van der Leun, A. M., Traets, J. J. H., Vos, J. L., Elbers, J. B. W., Patiwael, S., Qiao, X., Machuca-Ostos, M., Thommen, D. S., Haanen, J. B. A. G., Schumacher, T. N. M., & Zuur, C. L. (2023). Dual immune checkpoint blockade induces analogous alterations in the dysfunctional CD8+ T cell and activated Treg compartment. Cancer discovery, CD-22-0851. Advance online publication. https://doi.org/10.1158/2159-8290.CD-22-0851

Baran, Y., Bercovich, A., Sebe-Pedros, A., Lubling, Y., Giladi, A., Chomsky, E., Meir, Z., Hoichman, M., Lifshitz, A., & Tanay, A. (2019). MetaCell: analysis of single-cell RNA-seq data using K-nn graph partitions. Genome biology, 20(1), 206. https://doi.org/10.1186/s13059-019-1812-2

Kürten, C. H. L., Kulkarni, A., Cillo, A. R., Santos, P. M., Roble, A. K., Onkar, S., Reeder, C., Lang, S., Chen, X., Duvvuri, U., Kim, S., Liu, A., Tabib, T., Lafyatis, R., Feng, J., Gao, S. J., Bruno, T. C., Vignali, D. A. A., Lu, X., Bao, R., … Ferris, R. L. (2021). Investigating immune and non-immune cell interactions in head and neck tumors by single-cell RNA sequencing. Nature communications, 12(1), 7338. https://doi.org/10.1038/s41467-021-27619-4


<!-- README.md is generated from README.Rmd. Please edit that file -->
# CIDER

<!-- badges: start -->
<!-- badges: end -->
Clustering Single-cell RNA-Seq (scRNA-Seq) data from multiple samples or conditions are often challenged by confounding factors, such as batch effects and biologically relevant variability. Existing batch effect removal methods typically require strong assumptions on the composition of cell populations being near identical across samples. Here we present **CIDER**, a **meta-clustering workflow** based on inter-group similarity measures.

CIDER can:

1.  address the **clustering** task for confounded scRNA-Seq data, or
2.  assess the biological correctness of integration as **a test metric**, while it does **not** require the existence of prior cellular annotations.

For more informtion please see our [preprint](https://doi.org/10.1101/2021.03.29.437525).

## Installation

You can install CIDER from [github](https://github.com/zhiyhu/CIDER/) with:

``` r
# install.packages("devtools")
devtools::install_github('zhiyhu/CIDER')
```

## CIDER as an evaluation metric - Quick start

If you have scRNA-Seq data corrected by an integration algorithm (e.g. Seurat-CCA, Harmony, Scanrama...). You can use CIDER to evaluate if the biological populations are correctly aligned.

Before running CIDER evaluation functions, make sure that you have a Seurat object (e.g. `seu.integrated`) with corrected PCs in `seu.integrated@reductions$pca@cell.embeddings`. Seurat-CCA automatically put the corrected PCs there. If other methods are used, the corrected PCs can be added using `seu.integrated@reductions$pca@cell.embeddings <- corrected.PCs`.

``` r
library(CIDER)
seu.integrated <- hdbscan.seurat(seu.integrated)
ider <- getIDEr(seu.integrated, verbose = FALSE)
seu.integrated <- estimateProb(seu.integrated, ider)
```

The evaluation scores (IDER-based similarity and empirical p values) can be visualised by the `scatterPlot` function. See [here](https://zhiyhu.github.io/CIDER/articles/evaluation.html) for the full tutorial of using CIDER evaluation functions

``` r
p1 <- scatterPlot(seu.integrated, "tsne", colour.by = "similarity")
p2 <- scatterPlot(seu.integrated, "tsne", colour.by = "pvalue") 
plot_grid(p1,p2, ncol = 3)
```

![](man/figures/evaluation_scatterplot.png)

## Use CIDER for clustering tasks

![](man/figures/clustering_diagram.png)

### Quick start - asCIDER

``` r
ider <- getIDEr(seu, 
                group.by.var = "initial_cluster",
                batch.by.var = "Batch")
seu <- finalClustering(seu, ider, cutree.h = 0.45)
```

### Quick start - dnCIDER

<!--- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>.--->
## Citation

Z. Hu, A. A. Ahmed, C. Yau. An interpretable meta-clustering framework for single-cell RNA-Seq data integration and evaluation. *bioRxiv* 2021.03.29.437525; doi: <https://doi.org/10.1101/2021.03.29.437525>

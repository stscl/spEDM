# spEDM

[![spEDM website:
https://stscl.github.io/spEDM/](reference/figures/spEDM.png)](https://stscl.github.io/spEDM/)

***Sp**atial **E**mpirical **D**ynamic **M**odeling*

`spEDM` is an R package for spatial causal discovery. It extends
Empirical Dynamic Modeling (EDM) from time series to spatial
cross-sectional data, supports both vector and raster spatial formats
via integration with `sf` and `terra`, and enables data-driven causal
inference from spatial snapshots.

> *Refer to the package documentation <https://stscl.github.io/spEDM/>
> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=spEDM) with:

``` r
install.packages("spEDM", dep = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/spEDM) with:

``` r
install.packages("spEDM",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dep = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/spEDM)
  with:

``` r
if (!requireNamespace("devtools")) {
    install.packages("devtools")
}
devtools::install_github("stscl/spEDM",
                         build_vignettes = TRUE,
                         dep = TRUE)
```

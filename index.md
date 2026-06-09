# spEDM

[![spEDM website:
https://stscl.github.io/spEDM/](reference/figures/spEDM.png)](https://stscl.github.io/spEDM/)

***Sp**atial **E**mpirical **D**ynamic **M**odeling*

*spEDM* is an R package for spatial causal discovery. It extends
Empirical Dynamic Modeling (EDM) from time series to spatial
cross-sectional data, provides seamless support for vector and raster
spatial data via tight integration with the
[*sf*](https://CRAN.R-project.org/package=sf) and
[*terra*](https://CRAN.R-project.org/package=terra) packages, and
enables data-driven causal inference from spatial snapshots.

> *Refer to the package documentation <https://stscl.github.io/spEDM/>
> for more detailed information.*

## Installation

- Install from [CRAN](https://CRAN.R-project.org/package=spEDM) with:

``` r

install.packages("spEDM", dependencies = TRUE)
```

- Install binary version from
  [R-universe](https://stscl.r-universe.dev/spEDM) with:

``` r

install.packages("spEDM",
                 repos = c("https://stscl.r-universe.dev",
                           "https://cloud.r-project.org"),
                 dependencies = TRUE)
```

- Install from source code on [GitHub](https://github.com/stscl/spEDM)
  with:

``` r

if (!requireNamespace("pak", quietly = TRUE)) {
    install.packages("pak")
}
pak::pak("stscl/spEDM", dependencies = TRUE)
```

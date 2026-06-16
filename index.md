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

## CITATION

Please cite **[spEDM](https://doi.org/10.1080/13658816.2026.2687121)**
as:

    Lyu, W., Dai, S., Song, Y., Zhao, W., Yi, W., Xiao, Y., Jia, N., 2026. Measuring causal strengths from spatial cross-sectional data with geographical cross mapping cardinality. International Journal of Geographical Information Science 1–23. https://doi.org/10.1080/13658816.2026.2687121

A BibTeX entry for LaTeX users is:

``` bib
@article{lyu2026gcmc, 
    title = {Measuring causal strengths from spatial cross-sectional data with geographical cross mapping cardinality}, 
    ISSN = {1362-3087}, 
    DOI = {10.1080/13658816.2026.2687121}, 
    journal = {International Journal of Geographical Information Science}, 
    publisher = {Informa UK Limited}, 
    author = {Lyu, Wenbo and Dai, Shaoqing and Song, Yongze and Zhao, Wufan and Yi, Wen and Xiao, Yumiao and Jia, Nan}, 
    year = {2026}, 
    month = {June}, 
    pages = {1–23} 
}
```

## Reference

Lyu, W., Dai, S., Song, Y., Zhao, W., Yi, W., Xiao, Y., Jia, N., 2026.
Measuring causal strengths from spatial cross-sectional data with
geographical cross mapping cardinality. International Journal of
Geographical Information Science 1–23.
<https://doi.org/10.1080/13658816.2026.2687121>.

Lyu, W., Lei, Y., Yi, W., Song, Y., Li, X., Dai, S., Qin, Y., Zhao, W.,
2026. Causal discovery in urban data with temporal empirical dynamic
modeling: The R package tEDM. Computers, Environment and Urban Systems
127, 102435. <https://doi.org/10.1016/j.compenvurbsys.2026.102435>.

Gao, B., Yang, J., Chen, Z., Sugihara, G., Li, M., Stein, A., Kwan,
M.-P., Wang, J., 2023. Causal inference from cross-sectional earth
system data with geographical convergent cross mapping. Nature
Communications 14. <https://doi.org/10.1038/s41467-023-41619-6>.

Herrera, M., Mur, J., Ruiz, M., 2016. Detecting causal relationships
between spatial processes. Papers in Regional Science 95, 577–595.
<https://doi.org/10.1111/pirs.12144>.

Sugihara, G., May, R., Ye, H., Hsieh, C., Deyle, E., Fogarty, M., Munch,
S., 2012. Detecting Causality in Complex Ecosystems. Science 338,
496–500. <https://doi.org/10.1126/science.1227079>.

 

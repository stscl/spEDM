# Geographical Pattern Causality

## Methodological Background

Geographical Pattern Causality (GPC) infers causal relations from
spatial cross-sectional data by reconstructing a symbolic approximation
of the underlying spatial dynamical system.

Let \\x(s)\\ and \\y(s)\\ denote two spatial cross-sections over spatial
units \\s \in \mathcal{S}\\.

**(1) Spatial Embedding**

For each spatial unit \\s_i\\, GPC constructs an embedding vector

\\ \mathbf{E}\_{x(s_i)} = \big( x(s\_{i}^{(1)}), x(s\_{i}^{(2)}), \dots,
x(s\_{i}^{(E\tau)}) \big), \\

where \\s\_{i}^{(k)}\\ denotes the \\k\\-th spatially lagged value of
the spatial unit \\s_i\\, determined by embedding dimension \\E\\ and
lag \\\tau\\. This yields two reconstructed state spaces
\\\mathcal{M}\_x, \mathcal{M}\_y \subset \mathbb{R}^E\\.

**(2) Symbolic Pattern Extraction**

Local geometric transitions in each manifold are mapped to symbols

\\ \sigma_x(s_i),; \sigma_y(s_i) \in \mathcal{A}, \\

encoding increasing, decreasing, or non-changing modes. These symbolic
trajectories summarize local pattern evolution.

**(3) Cross-Pattern Mapping**

Causality from \\x \to y\\ is assessed by predicting:

\\ \hat{\sigma}\_y(s_i) = F\big( \sigma_x(s_j): s_j \in
\mathcal{N}\_k(s_i) \big), \\

where \\\mathcal{N}\_k\\ denotes the set of \\k\\ nearest neighbors in
\\\mathcal{M}\_x\\. The agreement structure between
\\\hat{\sigma}\_y(s_i)\\ and \\\sigma_y(s_i)\\ determines the causal
mode:

- Positive: \\\hat{\sigma}\_y = \sigma_y\\
- Negative: \\\hat{\sigma}\_y = -\sigma_y\\
- Dark: neither agreement nor opposition

**(4) Causal Strength**

The global causal strength is the normalized consistency of symbol
matches:

\\ C\_{x \to y} = \frac{1}{\|\mathcal{S}\|} \sum\_{s_i \in \mathcal{S}}
\mathbb{I}\big\[ \hat{\sigma}\_y(s_i) \bowtie \sigma_y(s_i) \big\], \\

where \\\bowtie\\ encodes positive, negative, or dark matching rules.

## Usage examples

### Example of spatial vector data

Load the `spEDM` package and its columbus spatial analysis data:

``` r
library(spEDM)

columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
columbus
## Simple feature collection with 49 features and 6 fields
## Geometry type: POLYGON
## Dimension:     XY
## Bounding box:  xmin: 5.874907 ymin: 10.78863 xmax: 11.28742 ymax: 14.74245
## Projected CRS: Undefined Cartesian SRS with unknown unit
## # A tibble: 49 × 7
##    hoval   inc  crime  open plumb discbd                                    geom
##    <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>                               <POLYGON>
##  1  80.5 19.5  15.7   2.85  0.217   5.03 ((8.624129 14.23698, 8.5597 14.74245, …
##  2  44.6 21.2  18.8   5.30  0.321   4.27 ((8.25279 14.23694, 8.282758 14.22994,…
##  3  26.4 16.0  30.6   4.53  0.374   3.89 ((8.653305 14.00809, 8.81814 14.00205,…
##  4  33.2  4.48 32.4   0.394 1.19    3.7  ((8.459499 13.82035, 8.473408 13.83227…
##  5  23.2 11.3  50.7   0.406 0.625   2.83 ((8.685274 13.63952, 8.677577 13.72221…
##  6  28.8 16.0  26.1   0.563 0.254   3.78 ((9.401384 13.5504, 9.434411 13.69427,…
##  7  75    8.44  0.178 0     2.40    2.74 ((8.037741 13.60752, 8.062716 13.60452…
##  8  37.1 11.3  38.4   3.48  2.74    2.89 ((8.247527 13.58651, 8.2795 13.5965, 8…
##  9  52.6 17.6  30.5   0.527 0.891   3.17 ((9.333297 13.27242, 9.671007 13.27361…
## 10  96.4 13.6  34.0   1.55  0.558   4.33 ((10.08251 13.03377, 10.0925 13.05275,…
## # ℹ 39 more rows
```

The false nearest neighbours (FNN) method helps identify the appropriate
minimal embedding dimension for reconstructing the state space of a time
series or spatial cross-sectional data.

``` r
spEDM::fnn(columbus, "crime", E = 1:10, eps = stats::sd(columbus$crime))
##        E:1        E:2        E:3        E:4        E:5        E:6        E:7 
## 0.59183673 0.04081633 0.04081633 0.10204082 0.00000000 0.00000000 0.00000000 
##        E:8 
## 0.00000000
```

The false nearest neighbours (FNN) ratio decreased to approximately
0.001 when the embedding dimension E reached 7, and remained relatively
stable thereafter. Therefore, we adopted \\E = 7\\ as the minimal
embedding dimension for subsequent parameter search.

Then, search optimal parameters:

``` r
# determine the type of causality using correlation
stats::cor.test(columbus$hoval,columbus$crime)
## 
##  Pearson's product-moment correlation
## 
## data:  columbus$hoval and columbus$crime
## t = -4.8117, df = 47, p-value = 1.585e-05
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  -0.7366777 -0.3497978
## sample estimates:
##        cor 
## -0.5744867

# since the correlation is -0.574, negative causality is selected as the metric to maximize in the optimal parameter search
spEDM::pc(columbus, "hoval", "crime", E = 5:6, k = 7:10, tau = 1, maximize = "negative")
## The suggested E,k,tau for variable crime is 6, 9 and 1
```

Run geographical pattern causality analysis

``` r
spEDM::gpc(columbus, "hoval", "crime", E = 6, k = 9)
##       type  strength      direction
## 1 positive       NaN hoval -> crime
## 2 negative 0.1340069 hoval -> crime
## 3     dark 0.1043991 hoval -> crime
## 4 positive       NaN crime -> hoval
## 5 negative 0.6251773 crime -> hoval
## 6     dark 0.1468990 crime -> hoval
```

Convergence diagnostics

``` r
crime_convergence = spEDM::gpc(columbus, "hoval", "crime",
                               libsizes = seq(5, 45, by = 5),
                               E = 6, k = 9, progressbar = FALSE)
crime_convergence
##    libsizes     type   strength      direction
## 1        10 positive 0.00000000 hoval -> crime
## 2        15 positive 0.00000000 hoval -> crime
## 3        20 positive 0.00000000 hoval -> crime
## 4        25 positive 0.00000000 hoval -> crime
## 5        30 positive 0.00000000 hoval -> crime
## 6        35 positive 0.00000000 hoval -> crime
## 7        40 positive 0.00000000 hoval -> crime
## 8        45 positive 0.00000000 hoval -> crime
## 9        10 negative 0.00000000 hoval -> crime
## 10       15 negative 0.10275505 hoval -> crime
## 11       20 negative 0.10662265 hoval -> crime
## 12       25 negative 0.11589173 hoval -> crime
## 13       30 negative 0.13260719 hoval -> crime
## 14       35 negative 0.12811574 hoval -> crime
## 15       40 negative 0.13875013 hoval -> crime
## 16       45 negative 0.13400692 hoval -> crime
## 17       10     dark 0.01608876 hoval -> crime
## 18       15     dark 0.01861307 hoval -> crime
## 19       20     dark 0.02173713 hoval -> crime
## 20       25     dark 0.01813809 hoval -> crime
## 21       30     dark 0.04362907 hoval -> crime
## 22       35     dark 0.04623515 hoval -> crime
## 23       40     dark 0.06830303 hoval -> crime
## 24       45     dark 0.09094789 hoval -> crime
## 25       10 positive 0.00000000 crime -> hoval
## 26       15 positive 0.00000000 crime -> hoval
## 27       20 positive 0.00000000 crime -> hoval
## 28       25 positive 0.00000000 crime -> hoval
## 29       30 positive 0.00000000 crime -> hoval
## 30       35 positive 0.00000000 crime -> hoval
## 31       40 positive 0.00000000 crime -> hoval
## 32       45 positive 0.00000000 crime -> hoval
## 33       10 negative 0.00000000 crime -> hoval
## 34       15 negative 0.00000000 crime -> hoval
## 35       20 negative 0.00000000 crime -> hoval
## 36       25 negative 0.00000000 crime -> hoval
## 37       30 negative 0.00000000 crime -> hoval
## 38       35 negative 0.10376552 crime -> hoval
## 39       40 negative 0.30424357 crime -> hoval
## 40       45 negative 0.62517731 crime -> hoval
## 41       10     dark 0.07618520 crime -> hoval
## 42       15     dark 0.08637332 crime -> hoval
## 43       20     dark 0.09795915 crime -> hoval
## 44       25     dark 0.10459567 crime -> hoval
## 45       30     dark 0.11259679 crime -> hoval
## 46       35     dark 0.11467853 crime -> hoval
## 47       40     dark 0.13381611 crime -> hoval
## 48       45     dark 0.14525108 crime -> hoval
plot(crime_convergence, ylimits = c(-0.01,1),
     xlimits = c(9,46), xbreaks = seq(10, 45, 10))
```

![Figure 1. Convergence curves of causal strengths among house value and
crime.](../reference/figures/gpc/fig1-1.png)

**Figure 1**. **Convergence curves of causal strengths among house value
and crime.**

  

### Example of spatial raster data

Load the `spEDM` package and its farmland NPP data:

``` r
library(spEDM)

npp = terra::rast(system.file("case/npp.tif", package = "spEDM"))
# To save the computation time, we will aggregate the data by 3 times
npp = terra::aggregate(npp, fact = 3, na.rm = TRUE)
npp
## class       : SpatRaster 
## size        : 135, 161, 5  (nrow, ncol, nlyr)
## resolution  : 30000, 30000  (x, y)
## extent      : -2625763, 2204237, 1867078, 5917078  (xmin, xmax, ymin, ymax)
## coord. ref. : CGCS2000_Albers 
## source(s)   : memory
## names       :      npp,        pre,      tem,      elev,         hfp 
## min values  :   187.50,   390.3351, -47.8194, -110.1494,  0.04434316 
## max values  : 15381.89, 23734.5330, 262.8576, 5217.6431, 42.68803711

# Check the validated cell number
nnamat = terra::as.matrix(npp[[1]], wide = TRUE)
nnaindice = which(!is.na(nnamat), arr.ind = TRUE)
dim(nnaindice)
## [1] 6920    2
```

Determining minimal embedding dimension:

``` r
spEDM::fnn(npp, "npp", E = 1:15,
           eps = stats::sd(terra::values(npp[["npp"]]),na.rm = TRUE))
##          E:1          E:2          E:3          E:4          E:5          E:6 
## 0.9813070569 0.5309427415 0.1322254335 0.0167630058 0.0017341040 0.0000000000 
##          E:7          E:8          E:9         E:10         E:11         E:12 
## 0.0001445087 0.0000000000 0.0000000000 0.0000000000 0.0002890173 0.0000000000 
##         E:13         E:14 
## 0.0001445087 0.0001445087
```

At \\E = 6\\, the false nearest neighbor ratio stabilizes approximately
at 0.0001 and remains constant thereafter. Therefore, \\E = 6\\ is
selected as minimal embedding dimension for the subsequent GPC analysis.

Then, search optimal parameters:

``` r
stats::cor.test(~ pre + npp,
                data = terra::values(npp[[c("pre","npp")]],
                                      datafame = TRUE, na.rm = TRUE))
## 
##  Pearson's product-moment correlation
## 
## data:  pre and npp
## t = 108.74, df = 6912, p-value < 2.2e-16
## alternative hypothesis: true correlation is not equal to 0
## 95 percent confidence interval:
##  0.7855441 0.8029426
## sample estimates:
##       cor 
## 0.7944062


g1 = spEDM::pc(npp, "npp", "pre", E = 6:10, k = 12, tau = 1:5, maximize = "positive")
g1
## The suggested E,k,tau for variable pre is 8, 12 and 5
```

Run geographical pattern causality analysis

``` r
spEDM::gpc(npp, "pre", "npp", E = 8, k = 12, tau = 5)
##       type  strength  direction
## 1 positive 0.5348284 pre -> npp
## 2 negative       NaN pre -> npp
## 3     dark 0.4184010 pre -> npp
## 4 positive 0.4346310 npp -> pre
## 5 negative 0.0000000 npp -> pre
## 6     dark 0.3810809 npp -> pre
```

Convergence diagnostics

``` r
npp_convergence = spEDM::gpc(npp, "pre", "npp",
                             libsizes = matrix(rep(seq(10,80,10),2),ncol = 2),
                             E = 8, k = 12, tau = 5, progressbar = FALSE)
npp_convergence
##    libsizes     type   strength  direction
## 1       100 positive 0.12538196 pre -> npp
## 2       400 positive 0.18018497 pre -> npp
## 3       900 positive 0.24754501 pre -> npp
## 4      1600 positive 0.30981312 pre -> npp
## 5      2500 positive 0.37799443 pre -> npp
## 6      3600 positive 0.43029973 pre -> npp
## 7      4900 positive 0.46526482 pre -> npp
## 8      6400 positive 0.50607359 pre -> npp
## 9       100 negative 0.00000000 pre -> npp
## 10      400 negative 0.00000000 pre -> npp
## 11      900 negative 0.00000000 pre -> npp
## 12     1600 negative 0.00000000 pre -> npp
## 13     2500 negative 0.00000000 pre -> npp
## 14     3600 negative 0.00000000 pre -> npp
## 15     4900 negative 0.00000000 pre -> npp
## 16     6400 negative 0.00000000 pre -> npp
## 17      100     dark 0.06183032 pre -> npp
## 18      400     dark 0.11232502 pre -> npp
## 19      900     dark 0.17258621 pre -> npp
## 20     1600     dark 0.23518092 pre -> npp
## 21     2500     dark 0.29193486 pre -> npp
## 22     3600     dark 0.33645457 pre -> npp
## 23     4900     dark 0.37582609 pre -> npp
## 24     6400     dark 0.40879652 pre -> npp
## 25      100 positive 0.07934454 npp -> pre
## 26      400 positive 0.11817574 npp -> pre
## 27      900 positive 0.18403265 npp -> pre
## 28     1600 positive 0.26820265 npp -> pre
## 29     2500 positive 0.34026886 npp -> pre
## 30     3600 positive 0.38219004 npp -> pre
## 31     4900 positive 0.41126235 npp -> pre
## 32     6400 positive 0.42254648 npp -> pre
## 33      100 negative 0.00000000 npp -> pre
## 34      400 negative 0.00000000 npp -> pre
## 35      900 negative 0.00000000 npp -> pre
## 36     1600 negative 0.00000000 npp -> pre
## 37     2500 negative 0.00000000 npp -> pre
## 38     3600 negative 0.00000000 npp -> pre
## 39     4900 negative 0.00000000 npp -> pre
## 40     6400 negative 0.00000000 npp -> pre
## 41      100     dark 0.02919551 npp -> pre
## 42      400     dark 0.06567401 npp -> pre
## 43      900     dark 0.11937259 npp -> pre
## 44     1600     dark 0.18083090 npp -> pre
## 45     2500     dark 0.23808200 npp -> pre
## 46     3600     dark 0.28469842 npp -> pre
## 47     4900     dark 0.32666764 npp -> pre
## 48     6400     dark 0.36603708 npp -> pre
plot(npp_convergence, ylimits = c(-0.01,0.65),
     xlimits = c(0,6500), xbreaks = seq(100, 6400, 500))
```

![Figure 2. Convergence curves of causal strengths among precipitation
and NPP.](../reference/figures/gpc/fig2-1.png)

**Figure 2**. **Convergence curves of causal strengths among
precipitation and NPP.**

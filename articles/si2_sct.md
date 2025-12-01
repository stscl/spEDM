# Spatial Causality Test

## Methodological Background

Let us begin by defining two spatial cross-sectional series
\\\\x_s\\\_{s \in S}\\ and \\\\y_s\\\_{s \in S}\\, where \\S\\
represents the study area.

We use:

\\ \mathcal{X}\_W = \\ W_i x \mid W_i \in \mathcal{W}(x, y) \\ \\

\\ \mathcal{Y}\_W = \\ W_i y \mid W_i \in \mathcal{W}(x, y) \\ \\

to denote the sets of spatial lags of \\x\\ and \\y\\ given by all the
weighting matrices in \\\mathcal{W}(x, y)\\(that is the set of spatial
dependence structures between \\x\\ and \\y\\).

We say that \\\\x_s\\\_{s \in S}\\ does not cause \\\\y_s\\\_{s \in S}\\
**under the spatial structures** \\\mathcal{X}\_W\\ and
\\\mathcal{Y}\_W\\ if

\\ h\_{y \|\mathcal{Y}\_W}(m) = h\_{y \|
\mathcal{Y}\_W,\mathcal{X}\_W}(m) \\ A unilateral non-parametric test
can be applied to assess the spatial causality via the following null
hypothesis:

\\ H_0 : \\x_s\\\_{s \in S} \text{ does not cause } \\y_s\\\_{s \in S}
\text{ under the spatial structures } \mathcal{X}\_W \text{ and }
\mathcal{Y}\_W, \\

with the following statistic:

\\ \hat{\delta}(\mathcal{Y}\_W, \mathcal{X}\_W) = \hat{h}\_{y
\|\mathcal{Y}\_W}(m) - \hat{h}\_{y \| \mathcal{Y}\_W,\mathcal{X}\_W}(m)
\\

where \\\hat{h}\_\*(m)\\ is the estimated conditional symbolic entropy
using Shannon’s entropy with \\m-1\\ nearest neighbors. The alternative
is that the null hypothesis of is not true.

If \\\mathcal{X}\_W\\ does not contain extra information about \\y\\
then \\\hat{\delta}(\mathcal{Y}\_W, \mathcal{X}\_W) = 0\\, otherwise, \$
(\_W, \_W) \> 0\$.

\\h\_{y \|\mathcal{Y}\_W}(m)\\ measures the uncertainty of the
distribution of symbols of \\y\\, conditional to the symbols of its
spatial lag, \$ \_W\$. Moreover, \\h\_{y \|
\mathcal{Y}\_W,\mathcal{X}\_W}(m)\\ measures the uncertainty of the
distribution of symbols of \\y\\, conditional to the symbols of the
spatial lags of \\y\\, \\\mathcal{Y}\_W\\, and of \\x\\, \$ \_W\$. If
the second variable, \\x\\, indeed causes the first one then there
should be a significant decrease in the entropy, and the statistic
\\\hat{\delta}(\mathcal{Y}\_W, \mathcal{X}\_W)\\ will take on high
positive values. If there is only a spatial correlation, but not
causation, the difference between both entropies will be small. The
statistical significance of \\\hat{\delta}(\mathcal{Y}\_W,
\mathcal{X}\_W)\\ is assessed using spatial block bootstrap.

## Usage examples

### An example of spatial lattice data

Load the `spEDM` package and the columbus OH dataset:

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
##    hoval   inc  crime  open plumb discbd                                     geom
##    <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>                                <POLYGON>
##  1  80.5 19.5  15.7   2.85  0.217   5.03 ((8.624129 14.23698, 8.5597 14.74245, 8…
##  2  44.6 21.2  18.8   5.30  0.321   4.27 ((8.25279 14.23694, 8.282758 14.22994, …
##  3  26.4 16.0  30.6   4.53  0.374   3.89 ((8.653305 14.00809, 8.81814 14.00205, …
##  4  33.2  4.48 32.4   0.394 1.19    3.7  ((8.459499 13.82035, 8.473408 13.83227,…
##  5  23.2 11.3  50.7   0.406 0.625   2.83 ((8.685274 13.63952, 8.677577 13.72221,…
##  6  28.8 16.0  26.1   0.563 0.254   3.78 ((9.401384 13.5504, 9.434411 13.69427, …
##  7  75    8.44  0.178 0     2.40    2.74 ((8.037741 13.60752, 8.062716 13.60452,…
##  8  37.1 11.3  38.4   3.48  2.74    2.89 ((8.247527 13.58651, 8.2795 13.5965, 8.…
##  9  52.6 17.6  30.5   0.527 0.891   3.17 ((9.333297 13.27242, 9.671007 13.27361,…
## 10  96.4 13.6  34.0   1.55  0.558   4.33 ((10.08251 13.03377, 10.0925 13.05275, …
## # ℹ 39 more rows
```

Detect spatial causality among the variables *inc*, *crime*, and *hoval*
:

``` r
# house value and crime
sc.test(columbus, "hoval", "crime", k = 15)
## spatial causality test
## hoval -> crime: statistic = 1.114, p value = 0.549
## crime -> hoval: statistic = 1.555, p value = 0.008

# household income and crime
sc.test(columbus, "inc", "crime", k = 15)
## spatial causality test
## inc -> crime: statistic = 0.739, p value = 0.99
## crime -> inc: statistic = 0.946, p value = 0.762

# household income and house value
sc.test(columbus, "inc", "hoval", k = 15)
## spatial causality test
## inc -> hoval: statistic = 1.211, p value = 0.311
## hoval -> inc: statistic = 0.992, p value = 0.702
```

### An example of spatial grid data

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
```

Detect spatial causality among the variables *pre*, *tem*, and *npp* :

``` r
# precipitation and npp
sc.test(npp,"pre","npp",k = 30)
## spatial causality test
## pre -> npp: statistic = 0.227, p value = 0
## npp -> pre: statistic = 0.16, p value = 0

# temperature and npp
sc.test(npp,"tem","npp",k = 30)
## spatial causality test
## tem -> npp: statistic = 0.168, p value = 0
## npp -> tem: statistic = 0.169, p value = 0

# precipitation and temperature
sc.test(npp,"pre","tem",k = 30)
## spatial causality test
## pre -> tem: statistic = 0.19, p value = 0
## tem -> pre: statistic = 0.147, p value = 0
```

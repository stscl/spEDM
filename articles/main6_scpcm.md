# Spatially Convergent Partial Cross Mapping

## Model principles

The methodological details are pending peer review and will be made
available thereafter.

## Usage examples

### Example of spatial vector data

Load the `spEDM` package and columbus OH spatial analysis dataset:

``` r

if (!requireNamespace("spEDM")) install.packages("spEDM")
## Loading required namespace: spEDM

columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
columbus
## Simple feature collection with 49 features and 6 fields
## Geometry type: POLYGON
## Dimension:     XY
## Bounding box:  xmin: 5.874907 ymin: 10.78863 xmax: 11.28742 ymax: 14.74245
## Projected CRS: Undefined Cartesian SRS with unknown unit
## # A tibble: 49 × 7
##    hoval   inc  crime  open plumb discbd                        geom
##    <dbl> <dbl>  <dbl> <dbl> <dbl>  <dbl>                   <POLYGON>
##  1  80.5 19.5  15.7   2.85  0.217   5.03 ((8.624129 14.23698, 8.559…
##  2  44.6 21.2  18.8   5.30  0.321   4.27 ((8.25279 14.23694, 8.2827…
##  3  26.4 16.0  30.6   4.53  0.374   3.89 ((8.653305 14.00809, 8.818…
##  4  33.2  4.48 32.4   0.394 1.19    3.7  ((8.459499 13.82035, 8.473…
##  5  23.2 11.3  50.7   0.406 0.625   2.83 ((8.685274 13.63952, 8.677…
##  6  28.8 16.0  26.1   0.563 0.254   3.78 ((9.401384 13.5504, 9.4344…
##  7  75    8.44  0.178 0     2.40    2.74 ((8.037741 13.60752, 8.062…
##  8  37.1 11.3  38.4   3.48  2.74    2.89 ((8.247527 13.58651, 8.279…
##  9  52.6 17.6  30.5   0.527 0.891   3.17 ((9.333297 13.27242, 9.671…
## 10  96.4 13.6  34.0   1.55  0.558   4.33 ((10.08251 13.03377, 10.09…
## # ℹ 39 more rows
```

We demonstrate how spatial vector data can be used in SCPCM analysis
through a causal example examining the influences of the level of
burglary incidents in a neighbourhood on house values, with
neighbourhood household income included as a conditioning variable.

Determine minimum embedding dimensions:

``` r

spEDM::fnn(columbus,"crime",E = 1:10)
## [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
## [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
## [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
##        E:1        E:2        E:3        E:4        E:5        E:6 
## 0.79591837 0.53061224 0.63265306 0.51020408 0.12244898 0.04081633 
##        E:7        E:8 
## 0.00000000 0.00000000
spEDM::fnn(columbus,"hoval",E = 1:10)
## [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
## [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
## [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
##        E:1        E:2        E:3        E:4        E:5        E:6 
## 0.85714286 0.77551020 0.51020408 0.61224490 0.22448980 0.08163265 
##        E:7        E:8 
## 0.00000000 0.00000000
spEDM::fnn(columbus,"inc",E = 1:10)
## [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
## [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
## [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
##        E:1        E:2        E:3        E:4        E:5        E:6 
## 0.73469388 0.24489796 0.30612245 0.38775510 0.24489796 0.04081633 
##        E:7        E:8 
## 0.00000000 0.00000000
```

Self prediction for parameter turning:

``` r

spEDM::simplex(columbus,"crime","crime",E = 7:10,k=12)
## The suggested E,k,tau for variable crime is 8, 12 and 1
spEDM::simplex(columbus,"hoval","hoval",E = 7:10,k=12)
## The suggested E,k,tau for variable hoval is 7, 12 and 1
spEDM::simplex(columbus,"inc","inc",E = 7:10,k=12)
## The suggested E,k,tau for variable inc is 8, 12 and 1
```

Conduct SCPCM:

``` r

crime_hoval = spEDM::scpcm(data = columbus,
                           cause = "crime",
                           effect = "hoval",
                           conds = "inc",
                           libsizes = seq(5, 45, by = 5),
                           E = c(8,7,8),
                           k = 12,
                           progressbar = FALSE)
crime_hoval
## -------------------------------------- 
## ***partial cross mapping prediction*** 
## -------------------------------------- 
##   libsizes crime->hoval hoval->crime
## 1        5   0.08948688  0.022401607
## 2       10   0.14511673  0.062986079
## 3       15   0.16849249  0.080653365
## 4       20   0.20940427  0.040974217
## 5       25   0.23160147  0.026613757
## 6       30   0.24709326  0.001942294
## 7       35   0.29707812 -0.014879976
## 8       40   0.33591810 -0.028758126
## 9       45   0.36337200 -0.040204136
## 
## ------------------------------ 
## ***cross mapping prediction*** 
## ------------------------------ 
##   libsizes crime->hoval hoval->crime
## 1        5    0.2224861    0.1378071
## 2       10    0.3023880    0.1880328
## 3       15    0.4005255    0.2085278
## 4       20    0.4681692    0.2203727
## 5       25    0.5098130    0.2293076
## 6       30    0.5467173    0.2326005
## 7       35    0.5732757    0.2331852
## 8       40    0.5822545    0.2361204
## 9       45    0.5883504    0.2364234
```

Visualize the result:

``` r

if (!requireNamespace("cowplot")) install.packages("cowplot")
## Loading required namespace: cowplot

fig1a = plot(crime_hoval,partial = FALSE,ylimits = c(0.1,0.75))
fig1b = plot(crime_hoval,partial = TRUE,ylimits = c(-0.05,0.55))
fig1 = cowplot::plot_grid(fig1a,fig1b,ncol = 2,label_fontfamily = 'serif',
                          labels = paste0('(',letters[1:2],')'))
fig1
```

![Figure 1. The cross mapping between crime and house value. a
Crime–hoval causality without accounting for covariates. b Crime–hoval
causality controlling for household
income.](../reference/figures/scpcm/fig1-1.png)

**Figure 1**. The cross mapping between crime and house value. **a**
Crime–hoval causality without accounting for covariates. **b**
Crime–hoval causality controlling for household income.

  

### Example of spatial raster data

Load the `spEDM` package and simulate raster data with a cyclic
interaction structure \\x \rightarrow y \rightarrow z \rightarrow x\\:

``` r

if (!requireNamespace("fields")) install.packages("fields")
## Loading required namespace: fields
if (!requireNamespace("MASS")) install.packages("MASS")

sim_trispecies = \(nx,ny,seed = 123){
  grid = expand.grid(seq(0, 10, length.out = nx),
                     seq(0, 10, length.out = ny))
  cov.fun = \(d, range = 1.5, sill = 1) sill * exp(-d/range)
  dist.mat = fields::rdist(grid)
  cov.mat = cov.fun(dist.mat, range = 1.5, sill = 1)
  set.seed(seed)
  res = replicate(3, {
    MASS::mvrnorm(1, rep(0, nrow(grid)), cov.mat) |>
      #pmax(0) |>
      sdsfun::normalize_vector(0,1) |>
      matrix(nrow = nx, ncol = ny) |>
      terra::rast()
  }, simplify = FALSE)
  terra::rast(res)
}

species = sim_trispecies(20,20, seed = 42)
names(species) = c("x","y","z")

sim = spEDM::slm(species, x = "x", y = "y", z = "z", k = 4,
                 step = 15, transient = 1, threshold = Inf,
                 aggregate_fn = \(.x) .x[length(.x)],
                 alpha_x = 0.75, alpha_y = 0.78, alpha_z = 0.76,
                 beta_xy = 0.35, beta_xz = 0, beta_yx = 0, beta_yz = 0.35,
                 beta_zx = 0.35, beta_zy = 0)

terra::values(species[["x"]]) = sim$x
terra::values(species[["y"]]) = sim$y
terra::values(species[["z"]]) = sim$z
species
## class       : SpatRaster
## size        : 20, 20, 3  (nrow, ncol, nlyr)
## resolution  : 1, 1  (x, y)
## extent      : 0, 20, 0, 20  (xmin, xmax, ymin, ymax)
## coord. ref. : 
## source(s)   : memory
## names       :        x,        y,        z
## min values  : 0.591053, 0.599642, 0.525905
## max values  : 0.881712, 0.864149, 0.876286
```

Determine minimum embedding dimensions:

``` r

spEDM::fnn(species, "x")
## [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
## [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
## [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
##    E:1    E:2    E:3    E:4    E:5    E:6    E:7    E:8    E:9 
## 0.8225 0.1750 0.0125 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
spEDM::fnn(species, "y")
## [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
## [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
## [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
##    E:1    E:2    E:3    E:4    E:5    E:6    E:7    E:8    E:9 
## 0.8625 0.1900 0.0250 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
spEDM::fnn(species, "z")
## [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
## [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
## [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
##  E:1  E:2  E:3  E:4  E:5  E:6  E:7  E:8  E:9 
## 0.87 0.13 0.00 0.00 0.00 0.00 0.00 0.00 0.00
```

Self prediction for parameter turning:

``` r

s1 = spEDM::simplex(species, "x", "x", E = 4:10, k = 12, tau = 1)
s2 = spEDM::simplex(species, "y", "y", E = 4:10, k = 12, tau = 1)
s3 = spEDM::simplex(species, "z", "z", E = 4:10, k = 12, tau = 1)
list(s1,s2,s3)
## [[1]]
## The suggested E,k,tau for variable x is 4, 12 and 1 
## 
## [[2]]
## The suggested E,k,tau for variable y is 4, 12 and 1 
## 
## [[3]]
## The suggested E,k,tau for variable z is 4, 12 and 1
```

Investigate the causation between x and z, with y as control variables:

``` r

xz = spEDM::scpcm(species, "x", "z", "y", E = 4, k = 12,
                  libsizes = matrix(seq(50,400,50), ncol = 1),
                  progressbar = FALSE)
xz
## -------------------------------------- 
## ***partial cross mapping prediction*** 
## -------------------------------------- 
##   libsizes      x->z      z->x
## 1       50 0.4736296 0.4933618
## 2      100 0.5212682 0.5577015
## 3      150 0.5581261 0.6211235
## 4      200 0.5729440 0.6705097
## 5      250 0.5593462 0.6971944
## 6      300 0.5321784 0.7161755
## 7      350 0.5133112 0.7527040
## 8      400 0.5014350 0.7640829
## 
## ------------------------------ 
## ***cross mapping prediction*** 
## ------------------------------ 
##   libsizes      x->z      z->x
## 1       50 0.6138023 0.6273096
## 2      100 0.6427045 0.6735237
## 3      150 0.6758891 0.7162489
## 4      200 0.7067327 0.7559218
## 5      250 0.7282486 0.7846552
## 6      300 0.7402653 0.8008604
## 7      350 0.7495153 0.8101575
## 8      400 0.7554037 0.8147119
```

Visualize the result:

``` r

fig2a = plot(xz,partial = FALSE,ylimits = c(0.6,0.95))
fig2b = plot(xz,partial = TRUE,ylimits = c(0.45,0.95))
fig2 = cowplot::plot_grid(fig2a,fig2b,ncol = 2,label_fontfamily = 'serif',
                          labels = paste0('(',letters[1:2],')'))
fig2
```

![Figure 2. The cross mapping between x and z. a x–z causality without
accounting for y. b x–z causality controlling for
y.](../reference/figures/scpcm/fig2-1.png)

**Figure 2**. The cross mapping between x and z. **a** x–z causality
without accounting for y. **b** x–z causality controlling for y.

# multiview embedding forecast

multiview embedding forecast

## Usage

``` r
# S4 method for class 'sf'
multiview(
  data,
  column,
  target,
  nvar,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  dist.average = TRUE,
  top = NULL,
  threads = detectThreads(),
  detrend = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
multiview(
  data,
  column,
  target,
  nvar,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  dist.average = TRUE,
  top = NULL,
  threads = detectThreads(),
  detrend = TRUE,
  grid.coord = TRUE
)
```

## Arguments

- data:

  observation data.

- column:

  name of library variable.

- target:

  name of target variable.

- nvar:

  number of variable combinations.

- E:

  (optional) embedding dimensions.

- k:

  (optional) number of nearest neighbors used.

- tau:

  (optional) step of spatial lags.

- style:

  (optional) embedding style (`0` includes current state, `1` excludes
  it).

- stack:

  (optional) whether to stack embeddings.

- lib:

  (optional) libraries indices (input needed: `vector` - spatial vector,
  `matrix` - spatial raster).

- pred:

  (optional) predictions indices (input requirement same as `lib`).

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- dist.average:

  (optional) whether to average distance.

- top:

  (optional) number of reconstructions used in MVE forecast.

- threads:

  (optional) number of threads to use.

- detrend:

  (optional) whether to remove the linear trend.

- nb:

  (optional) neighbours list.

- grid.coord:

  (optional) whether to detrend using cell center coordinates (`TRUE`)
  or row/column numbers (`FALSE`).

## Value

A vector (when input is sf object) or matrix

## References

Ye H., and G. Sugihara, 2016. Information leverage in interconnected
ecosystems: Overcoming the curse of dimensionality. Science 353:922-925.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
# \donttest{
multiview(columbus,
          column = c("inc","crime","open","plumb","discbd"),
          target = "hoval", nvar = 3)
#>  [1]  18.8646054  22.4869617  21.1550738  -1.1253460  -9.1007114  -2.7903615
#>  [7]  -7.6201012  -6.6613986  -8.7202885  16.2907046 -13.4100870   5.5617166
#> [13]  -4.6242104  -6.0215133 -12.1983725  -9.3303243  14.4441917   2.2055448
#> [19]  -4.3641292  -6.0871013 -11.5193470   2.6235111  13.7585733  -8.7467063
#> [25]  -9.7543853  -4.1516025   3.8137697 -18.6397054 -12.1462767 -12.0658493
#> [31]  10.9218398  19.4799373  -8.1531852  16.5727910   0.7379578   8.5333685
#> [37] -15.5058094 -17.8177107  23.9793164  27.1593848  18.1874514  -2.9723211
#> [43]  -2.8697946  -8.5016198   1.2200571   5.7208429  15.1875343   6.5063233
#> [49]   0.2310141
# }
```

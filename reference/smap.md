# smap forecast

smap forecast

## Usage

``` r
# S4 method for class 'sf'
smap(
  data,
  column,
  target,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  dist.average = TRUE,
  theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3,
    4, 6, 8),
  threads = detectThreads(),
  detrend = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
smap(
  data,
  column,
  target,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  dist.average = TRUE,
  theta = c(0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 0.75, 1, 1.5, 2, 3,
    4, 6, 8),
  threads = detectThreads(),
  detrend = TRUE,
  grid.coord = TRUE,
  embed.direction = 0
)
```

## Arguments

- data:

  observation data.

- column:

  name of library variable.

- target:

  name of target variable.

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

- theta:

  (optional) weighting parameter for distances.

- threads:

  (optional) number of threads to use.

- detrend:

  (optional) whether to remove the linear trend.

- nb:

  (optional) neighbours list.

- grid.coord:

  (optional) whether to detrend using cell center coordinates (`TRUE`)
  or row/column numbers (`FALSE`).

- embed.direction:

  (optional) direction selector for embeddings (`0` returns all
  directions, `1-8` correspond to NW, N, NE, W, E, SW, S, SE).

## Value

A list

- `xmap`:

  forecast performance

- `varname`:

  name of target variable

- `method`:

  method of cross mapping

## References

Sugihara G. 1994. Nonlinear forecasting for the classification of
natural time series. Philosophical Transactions: Physical Sciences and
Engineering, 348 (1688):477-495.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
smap(columbus,"inc","crime",E = 5,k = 6)
#> Warning: All evaluation metrics are identical within tolerance; choosing theta == 1 if available, otherwise closest to 1.
#> The suggested theta for variable crime is 1 
# }
```

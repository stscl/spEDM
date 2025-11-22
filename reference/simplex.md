# optimal parameter search for simplex forecasting

optimal parameter search for simplex forecasting

## Usage

``` r
# S4 method for class 'sf'
simplex(
  data,
  column,
  target,
  E = 2:10,
  k = E + 2,
  tau = 1,
  style = 1,
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  dist.average = TRUE,
  threads = detectThreads(),
  detrend = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
simplex(
  data,
  column,
  target,
  E = 2:10,
  k = E + 2,
  tau = 1,
  style = 1,
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  dist.average = TRUE,
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

Sugihara G. and May R. 1990. Nonlinear forecasting as a way of
distinguishing chaos from measurement error in time series. Nature,
344:734-741.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
# \donttest{
simplex(columbus,"inc","crime")
#> The suggested E,k,tau for variable crime is 5, 6 and 1 
# }
```

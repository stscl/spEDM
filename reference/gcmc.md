# geographical cross mapping cardinality

geographical cross mapping cardinality

## Usage

``` r
# S4 method for class 'sf'
gcmc(
  data,
  cause,
  effect,
  libsizes = NULL,
  E = 3,
  k = pmin(E^2),
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  threads = detectThreads(),
  detrend = FALSE,
  parallel.level = "low",
  bidirectional = TRUE,
  progressbar = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
gcmc(
  data,
  cause,
  effect,
  libsizes = NULL,
  E = 3,
  k = pmin(E^2),
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  threads = detectThreads(),
  detrend = FALSE,
  parallel.level = "low",
  bidirectional = TRUE,
  progressbar = TRUE,
  grid.coord = TRUE
)
```

## Arguments

- data:

  observation data.

- cause:

  name of causal variable.

- effect:

  name of effect variable.

- libsizes:

  (optional) number of spatial units used (input needed: `vector` -
  spatial vector, `matrix` - spatial raster).

- E:

  (optional) embedding dimensions.

- k:

  (optional) number of nearest neighbors.

- tau:

  (optional) step of spatial lags.

- style:

  (optional) embedding style (`0` includes current state, `1` excludes
  it).

- lib:

  (optional) libraries indices (input requirement same as `libsizes`).

- pred:

  (optional) predictions indices (input requirement same as `libsizes`).

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- threads:

  (optional) number of threads to use.

- detrend:

  (optional) whether to remove the linear trend.

- parallel.level:

  (optional) level of parallelism, `low` or `high`.

- bidirectional:

  (optional) whether to examine bidirectional causality.

- progressbar:

  (optional) whether to show the progress bar.

- nb:

  (optional) neighbours list.

- grid.coord:

  (optional) whether to detrend using cell center coordinates (`TRUE`)
  or row/column numbers (`FALSE`).

## Value

A list

- `xmap`:

  cross mapping results

- `cs`:

  causal strength

- `varname`:

  names of causal and effect variable

- `bidirectional`:

  whether to examine bidirectional causality

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
g = gcmc(columbus,"hoval","crime",E = 7,k = 18)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
g
#>   neighbors hoval->crime crime->hoval
#> 1        18     0.183642    0.2901235
# }
```

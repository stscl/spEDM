# intersection cardinality

intersection cardinality

## Usage

``` r
# S4 method for class 'sf'
ic(
  data,
  column,
  target,
  E = 2:10,
  k = E + 2,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  threads = detectThreads(),
  detrend = FALSE,
  parallel.level = "low",
  nb = NULL
)

# S4 method for class 'SpatRaster'
ic(
  data,
  column,
  target,
  E = 2:10,
  k = E + 2,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  threads = detectThreads(),
  detrend = FALSE,
  parallel.level = "low",
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

- E:

  (optional) embedding dimensions.

- k:

  (optional) number of nearest neighbors used.

- tau:

  (optional) step of spatial lags.

- style:

  (optional) embedding style (`0` includes current state, `1` excludes
  it).

- lib:

  (optional) libraries indices (input needed: `vector` - spatial vector,
  `matrix` - spatial raster).

- pred:

  (optional) predictions indices (input requirement same as `lib`).

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- threads:

  (optional) number of threads to use.

- detrend:

  (optional) whether to remove the linear trend.

- parallel.level:

  (optional) level of parallelism, `low` or `high`.

- nb:

  (optional) neighbours list.

- grid.coord:

  (optional) whether to detrend using cell center coordinates (`TRUE`)
  or row/column numbers (`FALSE`).

## Value

A list

- `xmap`:

  cross mapping performance

- `varname`:

  name of target variable

- `method`:

  method of cross mapping

- `tau`:

  step of time lag

## References

Tao, P., Wang, Q., Shi, J., Hao, X., Liu, X., Min, B., Zhang, Y., Li,
C., Cui, H., Chen, L., 2023. Detecting dynamical causality by
intersection cardinal concavity. Fundamental Research.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
ic(columbus,"hoval","crime", E = 7, k = 15:25)
#> The suggested E and k for variable crime is 7 and 18 
# }
```

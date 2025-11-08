# false nearest neighbours

false nearest neighbours

## Usage

``` r
# S4 method for class 'sf'
fnn(
  data,
  target,
  E = 1:10,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L1",
  rt = 10,
  eps = 2,
  threads = detectThreads(),
  detrend = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
fnn(
  data,
  target,
  E = 1:10,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L1",
  rt = 10,
  eps = 2,
  threads = detectThreads(),
  detrend = TRUE,
  grid.coord = TRUE
)
```

## Arguments

- data:

  observation data.

- target:

  name of target variable.

- E:

  (optional) embedding dimensions.

- tau:

  (optional) step of spatial lags.

- style:

  (optional) embedding style (`0` includes current state, `1` excludes
  it).

- lib:

  (optional) libraries indices.

- pred:

  (optional) predictions indices.

- dist.metric:

  (optional) distance metric (`L1`: Manhattan, `L2`: Euclidean).

- rt:

  (optional) escape factor.

- eps:

  (optional) neighborhood diameter.

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

A vector

## References

Kennel M. B., Brown R. and Abarbanel H. D. I., Determining embedding
dimension for phase-space reconstruction using a geometrical
construction, Phys. Rev. A, Volume 45, 3403 (1992).

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
fnn(columbus,"crime")
#>        E:1        E:2        E:3        E:4        E:5        E:6        E:7 
#> 0.79591837 0.53061224 0.63265306 0.51020408 0.12244898 0.04081633 0.00000000 
#>        E:8 
#> 0.00000000 
# }
```

# embedding spatial cross sectional data

embedding spatial cross sectional data

## Usage

``` r
# S4 method for class 'sf'
embedded(
  data,
  target,
  E = 3,
  tau = 1,
  style = 1,
  stack = FALSE,
  detrend = FALSE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
embedded(
  data,
  target,
  E = 3,
  tau = 1,
  style = 1,
  stack = FALSE,
  detrend = FALSE,
  grid.coord = TRUE,
  embed.direction = 0
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

- stack:

  (optional) whether to stack embeddings.

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

A matrix (when `stack` is `FALSE`) or list.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
v = embedded(columbus,"crime")
v[1:5,]
#>          [,1]     [,2]     [,3]
#> [1,] 24.71427 41.55964 43.45135
#> [2,] 26.24684 44.57868 40.73529
#> [3,] 29.41175 43.45135 37.55683
#> [4,] 34.64648 37.95647 38.69064
#> [5,] 40.46533 34.43484 42.28334

npp = terra::rast(system.file("case/npp.tif", package="spEDM"))
r = embedded(npp,"npp")
r[which(!is.na(r),arr.ind = TRUE)[1:5],]
#>          [,1]     [,2]     [,3]
#> [1,] 3467.572 3572.004 2611.333
#> [2,] 3355.386 3361.382      NaN
#> [3,] 3107.368 3736.406 3329.188
#> [4,] 2927.267 3519.752 3652.413
#> [5,] 3089.903 3562.450 3664.286
```

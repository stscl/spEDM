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
  stack = FALSE,
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
  stack = FALSE,
  lib = NULL,
  pred = NULL,
  dist.metric = "L1",
  rt = 10,
  eps = 2,
  threads = detectThreads(),
  detrend = TRUE,
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

- lib:

  (optional) libraries indices (input needed: `vector` - spatial vector,
  `matrix` - spatial raster).

- pred:

  (optional) predictions indices (input requirement same as `lib`).

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

- embed.direction:

  (optional) direction selector for embeddings (`0` returns all
  directions, `1-8` correspond to NW, N, NE, W, E, SW, S, SE).

## Value

A vector

## References

Kennel, M.B., Brown, R., Abarbanel, H.D.I., 1992. Determining embedding
dimension for phase-space reconstruction using a geometrical
construction. Physical Review A 45, 3403–3411.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
# \donttest{
spEDM::fnn(columbus,"crime")
#> [fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.
#> [fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.
#> [fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.
#>        E:1        E:2        E:3        E:4        E:5        E:6        E:7 
#> 0.79591837 0.53061224 0.63265306 0.51020408 0.12244898 0.04081633 0.00000000 
#>        E:8 
#> 0.00000000 
# }
```

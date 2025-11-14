# pattern causality

pattern causality

## Usage

``` r
# S4 method for class 'sf'
pc(
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
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  na.rm = TRUE,
  threads = detectThreads(),
  detrend = FALSE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
pc(
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
  zero.tolerance = max(k),
  relative = TRUE,
  weighted = TRUE,
  na.rm = TRUE,
  threads = detectThreads(),
  detrend = FALSE,
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

- zero.tolerance:

  (optional) maximum number of zeros tolerated in signature space.

- relative:

  (optional) whether to calculate relative changes in embeddings.

- weighted:

  (optional) whether to weight causal strength.

- na.rm:

  (optional) whether to remove `NA` samples in symbolic pattern
  generation.

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

A list

- `xmap`:

  cross mapping performance

- `varname`:

  name of target variable

- `method`:

  method of cross mapping

## References

Stavroglou, S.K., Pantelous, A.A., Stanley, H.E., Zuev, K.M., 2020.
Unveiling causal interactions in complex systems. Proceedings of the
National Academy of Sciences 117, 7599–7605.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
pc(columbus,"hoval","crime", E = 6:8, k = 7:10)
#> The suggested E,k,tau for variable crime is 8, 7 and 1 
# }
```

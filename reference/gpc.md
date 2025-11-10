# geographical pattern causality

geographical pattern causality

## Usage

``` r
# S4 method for class 'sf'
gpc(
  data,
  cause,
  effect,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  zero.tolerance = E + 2,
  relative = TRUE,
  weighted = TRUE,
  na.rm = TRUE,
  threads = detectThreads(),
  detrend = FALSE,
  bidirectional = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
gpc(
  data,
  cause,
  effect,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  dist.metric = "L2",
  zero.tolerance = E + 2,
  relative = TRUE,
  weighted = TRUE,
  na.rm = TRUE,
  threads = detectThreads(),
  detrend = FALSE,
  bidirectional = TRUE,
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

- E:

  (optional) embedding dimension.

- k:

  (optional) number of nearest neighbors used in cross mapping.

- tau:

  (optional) step of spatial lag.

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

- bidirectional:

  (optional) whether to examine bidirectional causality.

- nb:

  (optional) neighbours list.

- grid.coord:

  (optional) whether to detrend using cell center coordinates (`TRUE`)
  or row/column numbers (`FALSE`).

## Value

A list

- `causality`:

  per-sample causality statistics

- `summary`:

  overall causal strength

- `pattern`:

  pairwise pattern relationships

- `varname`:

  names of causal and effect variable

## References

Zhang, Z., Wang, J., 2025. A model to identify causality for geographic
patterns. International Journal of Geographical Information Science
1–21.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg", package="spEDM"))
# \donttest{
gpc(columbus,"hoval","crime", E = 3)
#> -------------------------------- 
#> ***pattern causality analysis*** 
#> -------------------------------- 
#>       type  strength      direction
#> 1 positive       NaN hoval -> crime
#> 2 negative 0.4672632 hoval -> crime
#> 3     dark 0.4354365 hoval -> crime
#> 4 positive 0.0000000 crime -> hoval
#> 5 negative 0.3920396 crime -> hoval
#> 6     dark 0.4048886 crime -> hoval
# }
```

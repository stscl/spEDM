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
#> $causality
#>    no positive  negative       dark     type valid
#> 1   0        0 0.0000000 0.37422083     dark  TRUE
#> 2   0        0 0.5129642 0.00000000 negative  TRUE
#> 3   0        0 0.6286459 0.00000000 negative  TRUE
#> 4   1        0 0.0000000 0.00000000       no  TRUE
#> 5   0        0 0.8673032 0.00000000 negative  TRUE
#> 6   0        0 0.0000000 0.99999999     dark  TRUE
#> 7   0        0 0.5186639 0.00000000 negative  TRUE
#> 8   1        0 0.0000000 0.00000000       no  TRUE
#> 9   0        0 0.3736861 0.00000000 negative  TRUE
#> 10  0        0 0.0000000 0.45286434     dark  TRUE
#> 11  0        0 0.0000000 0.99437287     dark  TRUE
#> 12  0        0 0.0000000 0.96865950     dark  TRUE
#> 13  0        0 0.6692985 0.00000000 negative  TRUE
#> 14  1        0 0.0000000 0.00000000       no  TRUE
#> 15  0        0 0.9543786 0.00000000 negative  TRUE
#> 16  0        0 0.9242713 0.00000000 negative  TRUE
#> 17  0        0 0.5051428 0.00000000 negative  TRUE
#> 18  0        0 0.0000000 0.61747155     dark  TRUE
#> 19  0        0 0.0000000 0.99948625     dark  TRUE
#> 20  0        0 0.0000000 0.72278605     dark  TRUE
#> 21  0        0 0.8138213 0.00000000 negative  TRUE
#> 22  1        0 0.0000000 0.00000000       no  TRUE
#> 23  0        0 0.5237650 0.00000000 negative  TRUE
#> 24  1        0 0.0000000 0.00000000       no  TRUE
#> 25  0        0 0.6723260 0.00000000 negative  TRUE
#> 26  1        0 0.0000000 0.00000000       no  TRUE
#> 27  0        0 0.0000000 0.64424426     dark  TRUE
#> 28  1        0 0.0000000 0.00000000       no  TRUE
#> 29  0        0 0.9569679 0.00000000 negative  TRUE
#> 30  1        0 0.0000000 0.00000000       no  TRUE
#> 31  0        0 0.0000000 0.13637658     dark  TRUE
#> 32  0        0 0.2732910 0.00000000 negative  TRUE
#> 33  1        0 0.0000000 0.00000000       no  TRUE
#> 34  0        0 0.0000000 0.08019356     dark  TRUE
#> 35  0        0 0.0000000 0.44506935     dark  TRUE
#> 36  1        0 0.0000000 0.00000000       no  TRUE
#> 37  0        0 0.9842995 0.00000000 negative  TRUE
#> 38  0        0 0.8381526 0.00000000 negative  TRUE
#> 39  0        0 0.3895475 0.00000000 negative  TRUE
#> 40  0        0 0.4364243 0.00000000 negative  TRUE
#> 41  1        0 0.0000000 0.00000000       no  TRUE
#> 42  0        0 0.0000000 0.12819493     dark  TRUE
#> 43  0        0 0.0000000 0.99235360     dark  TRUE
#> 44  1        0 0.0000000 0.00000000       no  TRUE
#> 45  0        0 0.6172039 0.00000000 negative  TRUE
#> 46  0        0 0.3642121 0.00000000 negative  TRUE
#> 47  0        0 0.3922613 0.00000000 negative  TRUE
#> 48  1        0 0.0000000 0.00000000       no  TRUE
#> 49  1        0 0.0000000 0.00000000       no  TRUE
#> 
#> $summary
#>       type  strength
#> 1 positive       NaN
#> 2 negative 0.4672632
#> 3     dark 0.4354365
#> 
#> $pattern
#>           11        13        31        33
#> 11       NaN 0.8087358       NaN 0.5330396
#> 13 0.0000000       NaN 0.5603749 0.3314576
#> 31 0.6560486 0.4150578       NaN 0.4961768
#> 33 0.3605805 0.5166237 0.2390131       NaN
#> 
# }
```

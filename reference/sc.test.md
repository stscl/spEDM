# spatial causality test

spatial causality test

## Usage

``` r
# S4 method for class 'sf'
sc.test(
  data,
  cause,
  effect,
  k,
  block = 3,
  boot = 399,
  seed = 42L,
  base = 2,
  lib = NULL,
  pred = NULL,
  threads = detectThreads(),
  detrend = TRUE,
  normalize = FALSE,
  progressbar = FALSE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
sc.test(
  data,
  cause,
  effect,
  k,
  block = 3,
  boot = 399,
  seed = 42L,
  base = 2,
  lib = NULL,
  pred = NULL,
  threads = detectThreads(),
  detrend = TRUE,
  normalize = FALSE,
  progressbar = FALSE,
  grid.coord = TRUE,
  embed.direction = 0
)
```

## Arguments

- data:

  observation data.

- cause:

  name of causal variable.

- effect:

  name of effect variable.

- k:

  (optional) number of nearest neighbors used in symbolization.

- block:

  (optional) number of blocks used in spatial block bootstrap.

- boot:

  (optional) number of bootstraps to perform.

- seed:

  (optional) random seed.

- base:

  (optional) logarithm base.

- lib:

  (optional) libraries indices (input requirement same as `libsizes`).

- pred:

  (optional) predictions indices (input requirement same as `libsizes`).

- threads:

  (optional) number of threads to use.

- detrend:

  (optional) whether to remove the linear trend.

- normalize:

  (optional) whether to normalize the result.

- progressbar:

  (optional) whether to show the progress bar.

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

- `sc`:

  statistic for spatial causality

- `varname`:

  names of causal and effect variables

## References

Herrera, M., Mur, J., Ruiz, M., 2016. Detecting causal relationships
between spatial processes. Papers in Regional Science 95, 577–595.

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
# \donttest{
spEDM::sc.test(columbus,"hoval","crime",k = 15)
#> spatial causality test
#> hoval -> crime: statistic = 1.114, p value = 0.549
#> crime -> hoval: statistic = 1.555, p value = 0.008
# }
```

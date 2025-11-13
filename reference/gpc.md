# geographical pattern causality

geographical pattern causality

## Usage

``` r
# S4 method for class 'sf'
gpc(
  data,
  cause,
  effect,
  libsizes = NULL,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  boot = 99,
  random = TRUE,
  seed = 42,
  dist.metric = "L2",
  zero.tolerance = k,
  relative = TRUE,
  weighted = TRUE,
  na.rm = TRUE,
  threads = detectThreads(),
  detrend = FALSE,
  parallel.level = "low",
  bidirectional = TRUE,
  progressbar = TRUE,
  nb = NULL
)

# S4 method for class 'SpatRaster'
gpc(
  data,
  cause,
  effect,
  libsizes = NULL,
  E = 3,
  k = E + 2,
  tau = 1,
  style = 1,
  lib = NULL,
  pred = NULL,
  boot = 99,
  random = TRUE,
  seed = 42,
  dist.metric = "L2",
  zero.tolerance = k,
  relative = TRUE,
  weighted = TRUE,
  na.rm = TRUE,
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

- boot:

  (optional) number of bootstraps to perform.

- random:

  (optional) whether to use random sampling.

- seed:

  (optional) random seed.

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

  cross mapping results (only present if `libsizes` is not `NULL`)

- `causality`:

  per-sample causality statistics(present if `libsizes` is `NULL`)

- `summary`:

  overall causal strength(present if `libsizes` is `NULL`)

- `pattern`:

  pairwise pattern relationships(present if `libsizes` is `NULL`)

- `varname`:

  names of causal and effect variable

- `bidirectional`:

  whether to examine bidirectional causality

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

# convergence diagnostics
g = gpc(columbus,"hoval","crime",libsizes = seq(5,45,5),E = 3)
#> Computing: [========================================] 100% (done)                         
#> Computing: [========================================] 100% (done)                         
g
#> -------------------------------- 
#> ***pattern causality analysis*** 
#> -------------------------------- 
#>    libsizes     type        mean        q05        q50       q95      direction
#> 1         5 positive 0.036176018 0.00000000 0.00000000 0.3575097 hoval -> crime
#> 2        10 positive 0.204854508 0.00000000 0.00000000 0.9998846 hoval -> crime
#> 3        15 positive 0.166966329 0.00000000 0.00000000 0.9745458 hoval -> crime
#> 4        20 positive 0.269291264 0.00000000 0.00000000 0.9900096 hoval -> crime
#> 5        25 positive 0.389912291 0.00000000 0.19954612 0.9990174 hoval -> crime
#> 6        30 positive 0.288498632 0.00000000 0.00000000 0.9977306 hoval -> crime
#> 7        35 positive 0.245224832 0.00000000 0.00000000 0.9525968 hoval -> crime
#> 8        40 positive 0.260125581 0.00000000 0.00000000 0.9808552 hoval -> crime
#> 9        45 positive 0.130895771 0.00000000 0.00000000 0.9410527 hoval -> crime
#> 10        5 negative 0.194028038 0.00000000 0.00000000 0.6936267 hoval -> crime
#> 11       10 negative 0.118665459 0.00000000 0.00000000 0.6023785 hoval -> crime
#> 12       15 negative 0.150688518 0.00000000 0.00000000 0.6119880 hoval -> crime
#> 13       20 negative 0.293656944 0.00000000 0.06997555 0.9712938 hoval -> crime
#> 14       25 negative 0.399200253 0.00000000 0.44550421 0.9803728 hoval -> crime
#> 15       30 negative 0.484154652 0.00000000 0.47930290 0.9763150 hoval -> crime
#> 16       35 negative 0.495635421 0.00000000 0.47393878 0.9831998 hoval -> crime
#> 17       40 negative 0.607135964 0.24102406 0.65585068 0.9859078 hoval -> crime
#> 18       45 negative 0.651939773 0.47901902 0.65604858 0.9379510 hoval -> crime
#> 19        5     dark 0.261691951 0.03075707 0.22042271 0.5848067 hoval -> crime
#> 20       10     dark 0.302543989 0.10982579 0.28860403 0.5315738 hoval -> crime
#> 21       15     dark 0.347725587 0.17785968 0.35142365 0.5243124 hoval -> crime
#> 22       20     dark 0.367281996 0.21419343 0.34947444 0.5273126 hoval -> crime
#> 23       25     dark 0.395459887 0.24837914 0.39863388 0.5263926 hoval -> crime
#> 24       30     dark 0.408233091 0.29256601 0.40037713 0.5344488 hoval -> crime
#> 25       35     dark 0.421067842 0.30944049 0.42995361 0.5190453 hoval -> crime
#> 26       40     dark 0.434209601 0.34181795 0.43099400 0.5305141 hoval -> crime
#> 27       45     dark 0.443724964 0.38512228 0.44136141 0.5095774 hoval -> crime
#> 28        5 positive 0.074643683 0.00000000 0.00000000 0.5622389 crime -> hoval
#> 29       10 positive 0.054697987 0.00000000 0.00000000 0.4455329 crime -> hoval
#> 30       15 positive 0.057658931 0.00000000 0.00000000 0.3431772 crime -> hoval
#> 31       20 positive 0.064427158 0.00000000 0.00000000 0.3871073 crime -> hoval
#> 32       25 positive 0.034785988 0.00000000 0.00000000 0.2948366 crime -> hoval
#> 33       30 positive 0.072017593 0.00000000 0.00000000 0.5058090 crime -> hoval
#> 34       35 positive 0.030759130 0.00000000 0.00000000 0.2666784 crime -> hoval
#> 35       40 positive 0.003632321 0.00000000 0.00000000 0.0000000 crime -> hoval
#> 36       45 positive 0.000000000 0.00000000 0.00000000 0.0000000 crime -> hoval
#> 37        5 negative 0.428359219 0.00000000 0.38208874 0.9553778 crime -> hoval
#> 38       10 negative 0.441918890 0.00000000 0.40488679 0.9684740 crime -> hoval
#> 39       15 negative 0.508595932 0.15587241 0.49527142 0.9030983 crime -> hoval
#> 40       20 negative 0.493236550 0.26067623 0.47004501 0.9200152 crime -> hoval
#> 41       25 negative 0.482211946 0.23931756 0.47234027 0.7902676 crime -> hoval
#> 42       30 negative 0.484994993 0.25158791 0.47416142 0.7412896 crime -> hoval
#> 43       35 negative 0.409372492 0.24545935 0.38740890 0.6320455 crime -> hoval
#> 44       40 negative 0.386542152 0.26427927 0.37003080 0.5630178 crime -> hoval
#> 45       45 negative 0.338634335 0.27048704 0.31378757 0.5019020 crime -> hoval
#> 46        5     dark 0.308907408 0.00000000 0.24524890 0.8626739 crime -> hoval
#> 47       10     dark 0.320420364 0.08725548 0.32736713 0.6176798 crime -> hoval
#> 48       15     dark 0.349842419 0.18894228 0.34680609 0.4907322 crime -> hoval
#> 49       20     dark 0.360706509 0.17112790 0.37015571 0.5115256 crime -> hoval
#> 50       25     dark 0.373620285 0.24108272 0.37437827 0.5366883 crime -> hoval
#> 51       30     dark 0.369014119 0.24340322 0.35727364 0.5041561 crime -> hoval
#> 52       35     dark 0.372023483 0.28005622 0.36791155 0.4661366 crime -> hoval
#> 53       40     dark 0.381208616 0.30244679 0.38098533 0.4670300 crime -> hoval
#> 54       45     dark 0.372072029 0.30867136 0.37282513 0.4425633 crime -> hoval
# }
```

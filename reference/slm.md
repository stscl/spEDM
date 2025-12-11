# spatial logistic map

spatial logistic map

## Usage

``` r
# S4 method for class 'sf'
slm(
  data,
  x = NULL,
  y = NULL,
  z = NULL,
  k = 4,
  step = 15,
  alpha_x = 0.28,
  alpha_y = 0.25,
  alpha_z = 0.22,
  beta_xy = 0.05,
  beta_xz = 0.05,
  beta_yx = 0.2,
  beta_yz = 0.2,
  beta_zx = 0.35,
  beta_zy = 0.35,
  threshold = Inf,
  transient = 1,
  interact = "local",
  aggregate_fn = NULL,
  noise = 0,
  seed = 42L,
  nb = NULL
)

# S4 method for class 'SpatRaster'
slm(
  data,
  x = NULL,
  y = NULL,
  z = NULL,
  k = 4,
  step = 15,
  alpha_x = 0.28,
  alpha_y = 0.25,
  alpha_z = 0.22,
  beta_xy = 0.05,
  beta_xz = 0.05,
  beta_yx = 0.2,
  beta_yz = 0.2,
  beta_zx = 0.35,
  beta_zy = 0.35,
  threshold = Inf,
  transient = 1,
  interact = "local",
  aggregate_fn = NULL,
  noise = 0,
  seed = 42L
)
```

## Arguments

- data:

  observation data.

- x:

  (optional) name of first spatial variable.

- y:

  (optional) name of second spatial variable.

- z:

  (optional) name of third spatial variable.

- k:

  (optional) number of neighbors to used.

- step:

  (optional) number of simulation time steps.

- alpha_x:

  (optional) growth parameter for x.

- alpha_y:

  (optional) growth parameter for y.

- alpha_z:

  (optional) growth parameter for z.

- beta_xy:

  (optional) cross-inhibition from x to y.

- beta_xz:

  (optional) cross-inhibition from x to z.

- beta_yx:

  (optional) cross-inhibition from y to x.

- beta_yz:

  (optional) cross-inhibition from y to z.

- beta_zx:

  (optional) cross-inhibition from z to x.

- beta_zy:

  (optional) cross-inhibition from z to y.

- threshold:

  (optional) set to `NaN` if the absolute value exceeds this threshold.

- transient:

  (optional) transients to be excluded from the results.

- interact:

  (optional) type of cross-variable interaction (`local` or
  `neighbors`).

- aggregate_fn:

  (optional) custom aggregation function (must accept a numeric vector
  and return a single numeric value).

- noise:

  (optional) standard deviation of white noise.

- seed:

  (optional) random seed.

- nb:

  (optional) neighbours list.

## Value

A list

## References

Willeboordse, F.H., The spatial logistic map as a simple prototype for
spatiotemporal chaos, Chaos, 533–540 (2003).

## Examples

``` r
columbus = sf::read_sf(system.file("case/columbus.gpkg",package="spEDM"))
columbus$inc = sdsfun::normalize_vector(columbus$inc)
slm(columbus,"inc")
#> $x
#>  [1] 0.8985452 0.8981079 0.8996913 0.9023158 0.9012219 0.8993983 0.9014626
#>  [8] 0.9011186 0.8995365 0.8998710 0.9016004 0.9012866 0.9013377 0.9012993
#> [15] 0.9013809 0.9016145 0.9011292 0.9007474 0.9009183 0.8921125 0.9012053
#> [22] 0.9011153 0.8965553 0.9003072 0.9015806 0.9015707 0.9009719 0.9016084
#> [29] 0.9015113 0.9003335 0.8988599 0.8961461 0.9011855 0.8993310 0.9005378
#> [36] 0.8983525 0.8996417 0.9010866 0.8981711 0.8921503 0.8963681 0.8953344
#> [43] 0.9005191 0.8994458 0.8999111 0.8982211 0.8981183 0.9004515 0.8985699
#> 
```

# spEDM 1.13

### enhancements

* Use queen contiguity to construct default `nb` when required but not provided (#1007).

* Integrate accepted IJGIS publication details for geographical cross mapping cardinality (#1005).

* No longer loads the `sf` and `terra` namespaces when loading `spEDM` package (#1000).

* Clarify output indexing and boundary handling in `fnn` generic via console output modification (#999).

### breaking changes

* Extend anisotropic spatial-lag embedding to fully support spatial raster data (#1012).

### bug fixes

* Fix incorrect auto-filling of default `E`, `tau`, `k` parameters in `scpcm` generic (#1009).

# spEDM 1.12

### new

* Correct maintainer surname spelling from `Lv` to `Lyu` for pinyin compliance (#966).

### enhancements

* Support data slicing by library and prediction sets for large scale pattern causality computation (#985).

* Assign zero to missing pattern causality values for positive, negative, and dark categories in returned results (#983).

* Qualify function calls in examples with namespace for reproducibility (#971).

* Inject plain-text terminal hint to `fnn` output (#965).

* Unify preprocessing strategies for embedding parameters `E` and `tau` across spatial causal discovery algorithms (#960).

* Enable asymmetric configurations (`E`, `tau`, `k`) in bidirectional pattern causality analysis (#955).

### breaking changes

* Match causality detection direction in optimal parameter search to actual variable naming (#995).

* Ensure that `k` is constrained to be greater than or equal to `E` during optimal parameter search (#990).

* Switch the `random` parameter to `replace` in the `gpc` generic implementation and revise the sampling logic accordingly (#987).

* Rename NA handling argument from `NA_rm` to `na_rm` and change default value to true (#968).

* Align `conds` parameter in `scpcm` generic to accept condition variables in direct path order for causal mediation analysis (#961).

### bug fixes

* Resolve visualization failure when probing only one direction of causality with `bidirectional` set to `FALSE` (#981).

* Correct mismatch between causal strength direction and prediction direction in pattern causality (#978).

* Fix incorrect midpoint-based pattern causality classification (#973).

# spEDM 1.11

### bug fixes

* Eliminate default initialization pollution from NaN-skipped elements in library and prediction sets (#951).

* Resolve macOS CRAN check failures caused by missing `#include` directives for `string` and `utility` headers (#948).

# spEDM 1.10

### new

* Provide R-level API and vignette for spatially convergent partial cross mapping (#936).

### enhancements

* Trim console output for geographical pattern causality results (#942).

* Enforce include order to prevent R headers from preceding Rcpp headers (#932).

* Support optional white noise in spatial logistic map simulation (#922).

* Align vignette metadata with index entry for correct cran listing (#920).

# spEDM 1.9

### new

* Provide R-level API and vignette for geographical pattern causality (#868).

* Introduce cross mapping of raster data with anisotropic embedding and variable sliding window steps (#812).

* Extend capabilities to incorporate mean-based and stacked spatial cross-sectional embeddings in cross mapping (#800).

### enhancements

* Drive extension of `slm` generics with custom aggregation capabilities (#912).

* Permit `simplex` and `ic` generics to accept varying E, k, and tau inputs (#900).

* Calibrate spacing conventions and enable real execution of available-thread detection (#885).

* Maintain consistent double precision floating-point comparisons across cpp sources (#850).

* Unify font specification in S3 plotting method for cross-mapping results (#830).

* Support spatial grid(raster) data detrending with cell center coordinates or row/column numbers (#815).

### breaking changes

* Remove redundant internal case data (`cu.tif`) from the package (#828).

* Harmonize parameter order across S4 generics for spatial vector and raster data (#817).

# spEDM 1.8

### new

* Provide new mechanisms for variable interaction in the spatial logistic map (#772).

* Introduce configurable distance metrics for cross mapping (#761).

* Enable alternative styles of spatial cross-sectional embedding (#755).

### enhancements

* Safeguard transient removal logic in spatial logistic map to prevent index errors (#744).

* Adjust the default `E` range in the `simplex` generic to `2:10` to support more robust reconstruction of state spaces (#739).

* Introduce multi-threading in distance-related computations where applicable to improve runtime efficiency (#718).

* Display p-value annotation for maximum library size in cross mapping visualization legend (#710).

### breaking changes

* Prevent unreliable predictions and potentially alter results when NaNs are present in cross mapping (#747).

* Rename coordinate columns in population density case csv to `lon` and `lat` (#741).

# spEDM 1.7

### new

* Provide R-level API and vignette for spatial logistic map ([#609](https://github.com/stscl/spEDM/pull/609),[#612](https://github.com/stscl/spEDM/pull/612)).

### enhancements

* No longer initializes the RNG state when loading `spEDM` package (#652).

* Replace logical vectors with integer index vectors for `lib` and `pred` in `simplex` and `s-mapping` forecasting sources (#632).

* Prevent duplicate generic registrations and simplify method definitions (#627).

* Unify parameter descriptions to lowercase (#570).

* Introduce confidence interval ribbon support in plot method for `gccm` results (#550).

### breaking changes

* Refine randomization strategy in spatial causality test (#643).

* Symbolization `C++` functions now compute medians from lib subset only (#599).

* Users must now use `detrend` instead of `trend.rm` (#559).

* Enable `column` parameter in `simplex()` and `smap()` and rename `columns` parameter to `column` in `multiview()` (#565).

### bug fixes

* Fix bug in plot method for `gccm` results where legend labels did not correctly match the corresponding line colors (#552).

# spEDM 1.6

### new

* Enact `fnn` R API support for false nearest neighbours method (#512).

* Integrate R API and vignette for geographical cross mapping cardinality (#455).

* Provide R-level API and vignette for spatial causality test (#403).

### enhancements

* Enable custom legend texts and colors in plot method for `gccm` results (#535).

* Reduce computational load in vignettes (#476).

* Document overall structure and usage of spEDM in a dedicated vignette (#415).

* Create `SSR` vignette for spatial cross-sectional data state-space reconstruction (#412).

* Include references for algorithms in `spEDM` (#367).

### breaking changes

* Use non-NA spatial units for `lib`/`pred` by default (#499).

* Refine internal case data (#348).

### bug fixes

* Patch memory error caused by mismatch between C++ (0-based) and R (1-based) indexing (#480).

* Fix error from non-matrix input in grid-type handling due to R matrix slicing (#474).

# spEDM 1.5

### new

* Enable `parallel.level` parameter to specify parallel granularity in `gccm` R API (#310).

* Implement the `multiview` function for multiview embedding forecasting (MVE) method (#221).

### enhancements

* Integrate `lib` parameter in `gccm` R API for library units selection (#278).

* Set the default `k` to `E+2` in the `gccm` R API (#261).

* Eliminate redundant computations at the source C++ code level (#233).

* Add `trend.rm` option in the R API for `embedded`, `simplex`, and `smap` methods to align with `gccm` behavior (#191).

* Refactor indexing of lag values and embedding vector generation for spatial lattice ([#186](https://github.com/stscl/spEDM/pull/186),[#184](https://github.com/stscl/spEDM/pull/184)) and grid data ([#183](https://github.com/stscl/spEDM/pull/183),[#181](https://github.com/stscl/spEDM/pull/181)).

### breaking changes

* Default plotting method places the legend in the top-left corner of the plot now (#325).

* Refine `simplex` & `smap` output on the R side (#263).

### bug fixes

* Fix bug in R functions `embedded`, `simplex`, `smap` when input data contains only one attribute column (#246).

# spEDM 1.4

### enhancements

* Improve default spatial neighbors list generation for spatial lattice data with support from the `sdsfun` package (#159).

### breaking changes

* Adjust the behavior of the `tau` parameter in the C++ source code and update the R side API (#154).

# spEDM 1.3

### new

* Implement the `smap` function to enable the selection of the optimal theta parameter (#128).

* Add `simplex` function to support selecting the optimal embedding dimension for variables (#98).

* Provide an R-level API for generating embeddings (#97).

### enhancements

* Now bidirectional mapping in the `gccm` result uses a `full join` structure when organized on the R side (#118).

* Support for calculating unidirectional mappings in the `gccm` function (#117).

* Relax `gccm` C++ source code `libsizes` minimum value constraint of `E+2` (#109).

* Provide a complete `GCCM` workflow for spatial lattice and grid data in the `gccm` vignette (#100).

* Support testing causal links in GCCM with different `E` and `k` for cause and effect variables (#96).

* Add thread settings for `gccm` generic (#94).

* Add `S-maps` cross-prediction support to `gccm` generic (#81).

### bug fixes

* Resolve r crash caused by invalid `E` [#90](https://github.com/stscl/spEDM/pull/90) and `k` [#89](https://github.com/stscl/spEDM/pull/89) parameter settings in `gccm`.

* Fix incorrect Pearson correlation calculation in `C++` code when input contains NA (#83).

# spEDM 1.2

### enhancements

* Encapsulate the `gccm` generic using the S4 class (#72).

* Add options for `tau`, `k`, and `progressbar` in `gccm` (#69).

* Add `print` and `plot` s3 methods for `gccm` result (#64).

### bug fixes

* Fix the bug where the `gccm` generic returns empty result when input grid data contains `NA` values (#61).

# spEDM 1.1

### bug fixes

* Resolve CRAN auto check issues, no significant API changes.

# spEDM 1.0

### new

* Implement `gccm` generic for spatial lattice and grid data using modern C++.

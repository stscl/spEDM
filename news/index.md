# Changelog

## spEDM 1.9

#### new

- Introduce cross mapping of raster data with anisotropic embedding and
  variable sliding window steps
  ([\#812](https://github.com/stscl/spEDM/issues/812)).

- Extend capabilities to incorporate mean-based and stacked spatial
  cross-sectional embeddings in cross mapping
  ([\#800](https://github.com/stscl/spEDM/issues/800)).

#### enhancements

- Maintain consistent double precision floating-point comparisons across
  cpp sources ([\#850](https://github.com/stscl/spEDM/issues/850)).

- Unify font specification in S3 plotting method for cross-mapping
  results ([\#830](https://github.com/stscl/spEDM/issues/830)).

- Support spatial grid(raster) data detrending with cell center
  coordinates or row/column numbers
  ([\#815](https://github.com/stscl/spEDM/issues/815)).

- Adjust `Makevars` to track latest rcpparmadillo release and armadillo
  updates improving build integration and forward compatibility
  ([\#785](https://github.com/stscl/spEDM/issues/785)).

#### breaking changes

- Remove redundant internal case data (`cu.tif`) from the package
  ([\#828](https://github.com/stscl/spEDM/issues/828)).

- Harmonize parameter order across s4 generics for spatial vector and
  raster data ([\#817](https://github.com/stscl/spEDM/issues/817)).

## spEDM 1.8

CRAN release: 2025-08-31

#### new

- Provide new mechanisms for variable interaction in the spatial
  logistic map ([\#772](https://github.com/stscl/spEDM/issues/772)).

- Introduce configurable distance metrics for cross mapping
  ([\#761](https://github.com/stscl/spEDM/issues/761)).

- Enable alternative styles of spatial cross-sectional embedding
  ([\#755](https://github.com/stscl/spEDM/issues/755)).

#### enhancements

- Safeguard transient removal logic in spatial logistic map to prevent
  index errors ([\#744](https://github.com/stscl/spEDM/issues/744)).

- Adjust the default `E` range in the `simplex` generic to `2:10` to
  support more robust reconstruction of state spaces
  ([\#739](https://github.com/stscl/spEDM/issues/739)).

- Introduce multi-threading in distance-related computations where
  applicable to improve runtime efficiency
  ([\#718](https://github.com/stscl/spEDM/issues/718)).

- Display p-value annotation for maximum library size in cross mapping
  visualization legend
  ([\#710](https://github.com/stscl/spEDM/issues/710)).

#### breaking changes

- Prevent unreliable predictions and potentially alter results when NaNs
  are present in cross mapping
  ([\#747](https://github.com/stscl/spEDM/issues/747)).

- Rename coordinate columns in population density case csv to `lon` and
  `lat` ([\#741](https://github.com/stscl/spEDM/issues/741)).

## spEDM 1.7

CRAN release: 2025-06-25

#### new

- Provide R-level API and vignette for spatial logistic map
  ([\#609](https://github.com/stscl/spEDM/pull/609),[\#612](https://github.com/stscl/spEDM/pull/612)).

#### enhancements

- Replace logical vectors with integer index vectors for `lib` and
  `pred` in `simplex` and `s-mapping` forecasting sources
  ([\#632](https://github.com/stscl/spEDM/issues/632)).

- Prevent duplicate generic registrations and simplify method
  definitions ([\#627](https://github.com/stscl/spEDM/issues/627)).

- Unify parameter descriptions to lowercase
  ([\#570](https://github.com/stscl/spEDM/issues/570)).

- Introduce confidence interval ribbon support in plot method for `gccm`
  results ([\#550](https://github.com/stscl/spEDM/issues/550)).

#### breaking changes

- Refine randomization strategy in spatial causality test
  ([\#643](https://github.com/stscl/spEDM/issues/643)).

- Symbolization `C++` functions now compute medians from lib subset only
  ([\#599](https://github.com/stscl/spEDM/issues/599)).

- Users must now use `detrend` instead of `trend.rm`
  ([\#559](https://github.com/stscl/spEDM/issues/559)).

- Enable `column` parameter in
  [`simplex()`](https://stscl.github.io/spEDM/reference/simplex.md) and
  [`smap()`](https://stscl.github.io/spEDM/reference/smap.md) and rename
  `columns` parameter to `column` in
  [`multiview()`](https://stscl.github.io/spEDM/reference/multiview.md)
  ([\#565](https://github.com/stscl/spEDM/issues/565)).

#### bug fixes

- Fix bug in plot method for `gccm` results where legend labels did not
  correctly match the corresponding line colors
  ([\#552](https://github.com/stscl/spEDM/issues/552)).

## spEDM 1.6

CRAN release: 2025-05-16

#### new

- Enact `fnn` R API support for false nearest neighbours method
  ([\#512](https://github.com/stscl/spEDM/issues/512)).

- Integrate R API and vignette for geographical cross mapping
  cardinality ([\#455](https://github.com/stscl/spEDM/issues/455)).

- Provide R-level API and vignette for spatial causality test
  ([\#403](https://github.com/stscl/spEDM/issues/403)).

#### enhancements

- Enable custom legend texts and colors in plot method for `gccm`
  results ([\#535](https://github.com/stscl/spEDM/issues/535)).

- Reduce computational load in vignettes
  ([\#476](https://github.com/stscl/spEDM/issues/476)).

- Document overall structure and usage of spEDM in a dedicated vignette
  ([\#415](https://github.com/stscl/spEDM/issues/415)).

- Create `SSR` vignette for spatial cross-sectional data state-space
  reconstruction ([\#412](https://github.com/stscl/spEDM/issues/412)).

- Include references for algorithms in `spEDM`
  ([\#367](https://github.com/stscl/spEDM/issues/367)).

#### breaking changes

- Use non-NA spatial units for `lib`/`pred` by default
  ([\#499](https://github.com/stscl/spEDM/issues/499)).

- Refine internal case data
  ([\#348](https://github.com/stscl/spEDM/issues/348)).

#### bug fixes

- Patch memory error caused by mismatch between C++ (0-based) and R
  (1-based) indexing
  ([\#480](https://github.com/stscl/spEDM/issues/480)).

- Fix error from non-matrix input in grid-type handling due to R matrix
  slicing ([\#474](https://github.com/stscl/spEDM/issues/474)).

## spEDM 1.5

CRAN release: 2025-03-12

#### new

- Enable `parallel.level` parameter to specify parallel granularity in
  `gccm` R API ([\#310](https://github.com/stscl/spEDM/issues/310)).

- Implement the `multiview` function for multiview embedding forecasting
  (MVE) method ([\#221](https://github.com/stscl/spEDM/issues/221)).

#### enhancements

- Integrate `lib` parameter in `gccm` R API for library units selection
  ([\#278](https://github.com/stscl/spEDM/issues/278)).

- Set the default `k` to `E+2` in the `gccm` R API
  ([\#261](https://github.com/stscl/spEDM/issues/261)).

- Eliminate redundant computations at the source C++ code level
  ([\#233](https://github.com/stscl/spEDM/issues/233)).

- Add `trend.rm` option in the R API for `embedded`, `simplex`, and
  `smap` methods to align with `gccm` behavior
  ([\#191](https://github.com/stscl/spEDM/issues/191)).

- Refactor indexing of lag values and embedding vector generation for
  spatial lattice
  ([\#186](https://github.com/stscl/spEDM/pull/186),[\#184](https://github.com/stscl/spEDM/pull/184))
  and grid data
  ([\#183](https://github.com/stscl/spEDM/pull/183),[\#181](https://github.com/stscl/spEDM/pull/181)).

#### breaking changes

- Default plotting method places the legend in the top-left corner of
  the plot now ([\#325](https://github.com/stscl/spEDM/issues/325)).

- Refine `simplex` & `smap` output on the R side
  ([\#263](https://github.com/stscl/spEDM/issues/263)).

#### bug fixes

- Fix bug in R functions `embedded`, `simplex`, `smap` when input data
  contains only one attribute column
  ([\#246](https://github.com/stscl/spEDM/issues/246)).

## spEDM 1.4

CRAN release: 2025-02-10

#### enhancements

- Improve default spatial neighbors list generation for spatial lattice
  data with support from the `sdsfun` package
  ([\#159](https://github.com/stscl/spEDM/issues/159)).

#### breaking changes

- Adjust the behavior of the `tau` parameter in the C++ source code and
  update the R side API
  ([\#154](https://github.com/stscl/spEDM/issues/154)).

## spEDM 1.3

CRAN release: 2025-02-01

#### new

- Implement the `smap` function to enable the selection of the optimal
  theta parameter ([\#128](https://github.com/stscl/spEDM/issues/128)).

- Add `simplex` function to support selecting the optimal embedding
  dimension for variables
  ([\#98](https://github.com/stscl/spEDM/issues/98)).

- Provide an R-level API for generating embeddings
  ([\#97](https://github.com/stscl/spEDM/issues/97)).

#### enhancements

- Now bidirectional mapping in the `gccm` result uses a `full join`
  structure when organized on the R side
  ([\#118](https://github.com/stscl/spEDM/issues/118)).

- Support for calculating unidirectional mappings in the `gccm` function
  ([\#117](https://github.com/stscl/spEDM/issues/117)).

- Relax `gccm` C++ source code `libsizes` minimum value constraint of
  `E+2` ([\#109](https://github.com/stscl/spEDM/issues/109)).

- Provide a complete `GCCM` workflow for spatial lattice and grid data
  in the `gccm` vignette
  ([\#100](https://github.com/stscl/spEDM/issues/100)).

- Support testing causal links in GCCM with different `E` and `k` for
  cause and effect variables
  ([\#96](https://github.com/stscl/spEDM/issues/96)).

- Add thread settings for `gccm` generic
  ([\#94](https://github.com/stscl/spEDM/issues/94)).

- Add `S-maps` cross-prediction support to `gccm` generic
  ([\#81](https://github.com/stscl/spEDM/issues/81)).

#### bug fixes

- Resolve r crash caused by invalid `E`
  [\#90](https://github.com/stscl/spEDM/pull/90) and `k`
  [\#89](https://github.com/stscl/spEDM/pull/89) parameter settings in
  `gccm`.

- Fix incorrect Pearson correlation calculation in `C++` code when input
  contains NA ([\#83](https://github.com/stscl/spEDM/issues/83)).

## spEDM 1.2

CRAN release: 2025-01-14

#### enhancements

- Encapsulate the `gccm` generic using the S4 class
  ([\#72](https://github.com/stscl/spEDM/issues/72)).

- Add options for `tau`, `k`, and `progressbar` in `gccm`
  ([\#69](https://github.com/stscl/spEDM/issues/69)).

- Add `print` and `plot` s3 methods for `gccm` result
  ([\#64](https://github.com/stscl/spEDM/issues/64)).

#### bug fixes

- Fix the bug where the `gccm` generic returns empty result when input
  grid data contains `NA` values
  ([\#61](https://github.com/stscl/spEDM/issues/61)).

## spEDM 1.1

CRAN release: 2024-12-19

#### bug fixes

- Resolve CRAN auto check issues, no significant API changes.

## spEDM 1.0

CRAN release: 2024-12-16

#### new

- Implement `gccm` generic for spatial lattice and grid data using
  modern C++.

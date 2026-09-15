#ifndef GCCM4Grid_H
#define GCCM4Grid_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <utility>
#include <limits>
#include <map>
#include "CppStats.h"
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.

/**
 * Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size and pred indice.
 *
 * This function calculates the cross mapping between a predictor variable (xEmbedings) and a response variable (yPred)
 * over a 2D grid, using either Simplex Projection or S-Mapping.
 *
 * @param xEmbedings           A 2D matrix of the predictor variable's embeddings (spatial cross-section data).
 * @param yPred                A 1D vector of the response variable's values (spatial cross-section data).
 * @param lib_sizes            A vector of two integers, where the first element is the row-wise library size and the second element is the column-wise library size.
 * @param possible_lib_indices A boolean vector indicating which spatial units are valid for inclusion in the library.
 * @param pred_indices         A vector indicating which spatial units to be predicted.
 * @param totalRow             The total number of rows in the 2D grid.
 * @param totalCol             The total number of columns in the 2D grid.
 * @param b                    The number of nearest neighbors to use for prediction.
 * @param simplex              If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta                The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads              The number of threads to use for parallel processing.
 * @param parallel_level       Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * @param row_size_mark        If true, use the row-wise libsize to mark the libsize; if false, use col-wise libsize.
 * @param dist_metric          Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average         Whether to average distance by the number of valid vector components.
 * @param win_ratios           Scale the sliding window step relative to the matrix width/height to speed up state-space predictions.
 *
 * @return  A vector of pairs, where each pair contains the library size and the corresponding cross mapping result.
 */
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<int>& lib_sizes,
    const std::vector<bool>& possible_lib_indices,
    const std::vector<int>& pred_indices,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool row_size_mark,
    int dist_metric,
    bool dist_average,
    const std::vector<double>& win_ratios = {0,0}
);

// Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size and pred indice (composite embeddings version).
std::vector<std::pair<int, double>> GCCMSingle4Grid(
    const std::vector<std::vector<std::vector<double>>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<int>& lib_sizes,
    const std::vector<bool>& possible_lib_indices,
    const std::vector<int>& pred_indices,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool row_size_mark,
    int dist_metric,
    bool dist_average,
    const std::vector<double>& win_ratios = {0,0}
);

/**
 * Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size.
 *
 * This function follows the same library construction logic as GCCMSingle4Lattice, where libraries
 * are created by selecting consecutive indices from possible_lib_indices with possible wraparound.
 *
 * @param xEmbedings           State-space embeddings for the predictor variable (each row is a spatial vector)
 * @param yPred                Target spatial cross-section series
 * @param lib_size             Number of consecutive spatial units to include in each library
 * @param lib_indices          Integer vector indicating the indices of eligible spatial units for library construction
 * @param pred_indices         Integer vector indicating spatial units to predict
 * @param totalRow             Total rows in spatial grid
 * @param totalCol             Total columns in spatial grid
 * @param b                    Number of nearest neighbors for prediction
 * @param simplex              Use simplex projection if true, S-mapping if false
 * @param theta                Distance weighting parameter for S-mapping
 * @param threads              The number of threads to use for parallel processing
 * @param parallel_level       Level of parallel computing: 0 for `lower`, 1 for `higher`
 * @param dist_metric          Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average         Whether to average distance by the number of valid vector components.
 *
 * @return A vector of pairs, where each pair contains the library size and the corresponding cross mapping result.
 */
std::vector<std::pair<int, double>> GCCMSingle4GridOneDim(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    int dist_metric,
    bool dist_average
);

// Perform Grid-based Geographical Convergent Cross Mapping (GCCM) for a single library size (composite embeddings version).
std::vector<std::pair<int, double>> GCCMSingle4GridOneDim(
    const std::vector<std::vector<std::vector<double>>>& xEmbedings,
    const std::vector<double>& yPred,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    int totalRow,
    int totalCol,
    int b,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    int dist_metric,
    bool dist_average
);

/**
 * Perform Geographical Convergent Cross Mapping (GCCM) for spatial grid data.
 *
 * This function calculates the cross mapping between predictor variables (xMatrix) and response variables (yMatrix)
 * over a 2D grid, using either Simplex Projection or S-Mapping. It supports parallel processing and progress tracking.
 *
 * @param xMatrix        A 2D matrix of the predictor variable's values (spatial cross-section data).
 * @param yMatrix        A 2D matrix of the response variable's values (spatial cross-section data).
 * @param lib_sizes      A 2D vector where the first sub-vector contains row-wise library sizes and the second sub-vector contains column-wise library sizes.
 * @param lib            A vector of pairs representing the indices (row, column) of spatial units to be the library.
 * @param pred           A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param E              The number of dimensions for attractor reconstruction.
 * @param tau            The step of spatial lags for prediction.
 * @param b              The number of nearest neighbors to use for prediction.
 * @param simplex        If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta          The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads        The number of threads to use for parallel processing.
 * @param parallel_level Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * @param style          Embedding style selector (0: includes current state, 1: excludes it).
 * @param stack          Embedding arrangement selector (0: single - average lags, 1: composite - stack).  Default is 0 (average lags).
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average   Whether to average distance by the number of valid vector components.
 * @param single_sig     Whether to estimate significance and confidence intervals using a single rho value.
 * @param progressbar    If true, display a progress bar during computation.
 * @param dir            Direction selector for embeddings where 0 returns all directions for embeddings, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 * @param win_rations    Scale the sliding window step relative to the matrix width/height to speed up state-space predictions.
 *
 * @return A 2D vector where each row contains the library size, mean cross mapping result,
 *         significance, and confidence interval bounds.
 */
std::vector<std::vector<double>> GCCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<std::vector<int>>& lib_sizes,
    const std::vector<std::pair<int, int>>& lib,
    const std::vector<std::pair<int, int>>& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    int parallel_level,
    int style,
    int stack,
    int dist_metric,
    bool dist_average,
    bool single_sig,
    bool progressbar = false,
    const std::vector<int>& dir = {0},
    const std::vector<double>& win_ratios = {0,0}
);

/**
 * Perform Geographical Convergent Cross Mapping (GCCM) for spatial grid data.
 *
 * This function calculates the cross mapping between predictor variables (xMatrix) and response variables (yMatrix)
 * over a 2D grid, using either Simplex Projection or S-Mapping. It supports parallel processing and progress tracking.
 *
 * @param xMatrix        A 2D matrix of the predictor variable's values (spatial cross-section data).
 * @param yMatrix        A 2D matrix of the response variable's values (spatial cross-section data).
 * @param lib_sizes      Number of consecutive spatial units to include in each library.
 * @param lib            A vector of representing the indices of spatial units to be the library.
 * @param pred           A vector of representing the indices of spatial units to be predicted.
 * @param E              The number of dimensions for attractor reconstruction.
 * @param tau            The step of spatial lags for prediction.
 * @param b              The number of nearest neighbors to use for prediction.
 * @param simplex        If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta          The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads        The number of threads to use for parallel processing.
 * @param parallel_level Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * @param style          Embedding style selector (0: includes current state, 1: excludes it).
 * @param stack          Embedding arrangement selector (0: single - average lags, 1: composite - stack).  Default is 0 (average lags).
 * @param dist_metric    Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average   Whether to average distance by the number of valid vector components.
 * @param single_sig     Whether to estimate significance and confidence intervals using a single rho value.
 * @param progressbar    If true, display a progress bar during computation.
 * @param dir            Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * @return A 2D vector where each row contains the library size, mean cross mapping result,
 *         significance, and confidence interval bounds.
 */
std::vector<std::vector<double>> GCCM4GridOneDim(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    int E,
    int tau,
    int b,
    bool simplex,
    double theta,
    int threads,
    int parallel_level,
    int style,
    int stack,
    int dist_metric,
    bool dist_average,
    bool single_sig,
    bool progressbar = false,
    const std::vector<int>& dir = {0}
);

#endif // GCCM4Grid_H

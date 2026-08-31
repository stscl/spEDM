#ifndef SCPCM4Grid_H
#define SCPCM4Grid_H

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
#include "DataStruct.h"
// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.

/**
 * @brief Computes the partial correlation between a spatial cross-sectional series and its prediction
 *        using the Simplex Projection method, incorporating control variables in a grid-based spatial setting.
 *
 * This function reconstructs the state-space and applies Simplex Projection prediction while accounting
 * for control variables in a spatial grid data. The process can be cumulative or independent in incorporating
 * control variables.
 *
 * @param vectors: Reconstructed state-space where each row represents a separate vector/state.
 * @param target: Spatial cross-sectional series used as the prediction target.
 * @param controls: Cross-sectional data of control variables, stored row-wise.
 * @param lib_indices: Integer vector indicating which states to include when searching for neighbors.
 * @param pred_indices: Integer vector indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for Simplex Projection.
 * @param nrow: Number of rows in the input spatial grid data.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @param style: Embedding style selector (0: includes current state, 1: excludes it).
 * @param dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average: Whether to average distance by the number of valid vector components.
 * @param dir: Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * @return A vector of size 2 containing:
 *         - rho[0]: Pearson correlation between the target and its predicted values.
 *         - rho[1]: Partial correlation between the target and its predicted values, adjusting for control variables.
 */
std::vector<double> PartialSimplex4Grid(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    int nrow,
    bool cumulate = false,
    int style = 1,
    int dist_metric = 2,
    bool dist_average = true,
    const std::vector<int>& dir = {0}
);

/**
 * Computes the partial correlation between a spatial cross-sectional series and its prediction
 * using the Simplex Projection method, incorporating control variables in a grid-based spatial
 * setting (composite embeddings version).
 */
std::vector<double> PartialSimplex4Grid(
    const std::vector<std::vector<std::vector<double>>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    int nrow,
    bool cumulate = false,
    int style = 1,
    int dist_metric = 2,
    bool dist_average = true,
    const std::vector<int>& dir = {0}
);

/**
 * @brief Computes the partial correlation between a spatial cross-sectional series and its prediction
 *        using the S-Map method, incorporating control variables in a grid-based spatial setting.
 *
 * This function reconstructs the state-space and applies S-Map prediction while accounting for
 * control variables in a spatial grid data. The process can be cumulative or independent in incorporating
 * control variables.
 *
 * @param vectors: Reconstructed state-space where each row represents a separate vector/state.
 * @param target: Spatial cross-sectional series used as the prediction target.
 * @param controls: Cross-sectional data of control variables, stored row-wise.
 * @param lib_indices: Integer vector indicating which states to include when searching for neighbors.
 * @param pred_indices: Integer vector indicating which states to predict from.
 * @param conEs: Vector specifying the number of dimensions for attractor reconstruction with control variables.
 * @param taus: Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param num_neighbors: Vector specifying the numbers of neighbors to use for Simplex Projection.
 * @param nrow: Number of rows in the input spatial grid data.
 * @param theta: Weighting parameter for distances in the S-Map method.
 * @param cumulate: Boolean flag to determine whether to cumulate the partial correlations.
 * @param style: Embedding style selector (0: includes current state, 1: excludes it).
 * @param dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average: Whether to average distance by the number of valid vector components.
 * @param dir: Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * @return A vector of size 2 containing:
 *         - rho[0]: Pearson correlation between the target and its predicted values.
 *         - rho[1]: Partial correlation between the target and its predicted values, adjusting for control variables.
 */
std::vector<double> PartialSMap4Grid(
    const std::vector<std::vector<double>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    int nrow,
    double theta = 1.0,
    bool cumulate = false,
    int style = 1,
    int dist_metric = 2,
    bool dist_average = true,
    const std::vector<int>& dir = {0}
);

/**
 * Computes the partial correlation between a spatial cross-sectional series and its prediction
 * using the S-Map method, incorporating control variables in a grid-based spatial setting (composite embeddings version).
 */
std::vector<double> PartialSMap4Grid(
    const std::vector<std::vector<std::vector<double>>>& vectors,
    const std::vector<double>& target,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& num_neighbors,
    int nrow,
    double theta = 1.0,
    bool cumulate = false,
    int style = 1,
    int dist_metric = 2,
    bool dist_average = true,
    const std::vector<int>& dir = {0}
);

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size.
 *
 * This function calculates the partial cross mapping between a predictor variable (xEmbedings) and a response
 * variable (yPred) over a 2D grid, using either Simplex Projection or S-Mapping.
 *
 * @param xEmbedings           A 2D matrix of the predictor variable's embeddings (spatial cross-sectional data).
 * @param yPred                A 1D vector of the response variable's values (spatial cross-sectional data).
 * @param controls             A 2D matrix that stores the control variables.
 * @param lib_sizes            A vector of two integers, where the first element is the row-wise library size and the second element is the column-wise library size.
 * @param possible_lib_indices A boolean vector indicating which spatial units are valid for inclusion in the library.
 * @param pred_indices         A integer vector indicating which spatial units to be predicted.
 * @param conEs                Number of dimensions for the attractor reconstruction with control variables
 * @param taus:                Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
 * @param totalRow             The total number of rows in the 2D grid.
 * @param totalCol             The total number of columns in the 2D grid.
 * @param b                    The numbers of nearest neighbors to use for prediction.
 * @param simplex              If true, use Simplex Projection; if false, use S-Mapping.
 * @param theta                The distance weighting parameter for S-Mapping (ignored if simplex is true).
 * @param threads              The number of threads to use for parallel processing.
 * @param parallel_level       Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * @param cumulate             Whether to cumulate the partial correlations.
 * @param row_size_mark        If ture, use the row-wise libsize to mark the libsize; if false, use col-wise libsize.
 * @param style                Embedding style selector (0: includes current state, 1: excludes it).
 * @param dist_metric          Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param dist_average         Whether to average distance by the number of valid vector components.
 * @param dir                  Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 * @param win_ratios           Scale the sliding window step relative to the matrix width/height to speed up state-space predictions.
 *
 * @return  A vector contains the library size and the corresponding cross mapping and partial cross mapping result.
 */
std::vector<PartialCorRes> SCPCMSingle4Grid(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_sizes,
    const std::vector<bool>& possible_lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    int totalRow,
    int totalCol,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool cumulate,
    bool row_size_mark,
    int style,
    int dist_metric,
    bool dist_average,
    const std::vector<int>& dir = {0},
    const std::vector<double>& win_ratios = {0,0}
);

// Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size (composite embeddings version).
std::vector<PartialCorRes> SCPCMSingle4Grid(
    const std::vector<std::vector<std::vector<double>>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    const std::vector<int>& lib_sizes,
    const std::vector<bool>& possible_lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    int totalRow,
    int totalCol,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool cumulate,
    bool row_size_mark,
    int style,
    int dist_metric,
    bool dist_average,
    const std::vector<int>& dir = {0},
    const std::vector<double>& win_ratios = {0,0}
);

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size.
 *
 * This function follows the same library construction logic as SCPCMSingle4Lattice, where libraries
 * are created by selecting consecutive indices from possible_lib_indices with possible wraparound.
 *
 * @param xEmbedings           State-space embeddings for the predictor variable (each row is a spatial vector)
 * @param yPred                Target spatial cross-sectional series
 * @param controls             Control variables stored by row
 * @param lib_size             Number of consecutive spatial units to include in each library
 * @param lib_indices          Integer vector indicating the indices of eligible spatial units for library construction
 * @param pred_indices         Integer vector indicating spatial units to predict
 * @param conEs                Embedding dimensions for control variables
 * @param taus                 Spatial lag steps for control variable embeddings
 * @param b                    Number of nearest neighbors for prediction
 * @param totalRow             Total rows in spatial grid
 * @param totalCol             Total columns in spatial grid
 * @param simplex              Use simplex projection if true, S-mapping if false
 * @param theta                Distance weighting parameter for S-mapping
 * @param threads              The number of threads to use for parallel processing
 * @param parallel_level       Level of parallel computing: 0 for `lower`, 1 for `higher`
 * @param cumulate             Enable cumulative partial correlations
 * @param style                Embedding style selector (0: includes current state, 1: excludes it)
 * @param dist_metric          Distance metric selector (1: Manhattan, 2: Euclidean)
 * @param dist_average         Whether to average distance by the number of valid vector components
 * @param dir                  Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * @return Vector of PartialCorRes containing mapping results for each library configuration
 */
std::vector<PartialCorRes> SCPCMSingle4GridOneDim(
    const std::vector<std::vector<double>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    int totalRow,
    int totalCol,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool cumulate,
    int style,
    int dist_metric,
    bool dist_average,
    const std::vector<int>& dir = {0}
);

// Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for a single library size (composite embeddings version).
std::vector<PartialCorRes> SCPCMSingle4GridOneDim(
    const std::vector<std::vector<std::vector<double>>>& xEmbedings,
    const std::vector<double>& yPred,
    const std::vector<std::vector<double>>& controls,
    int lib_size,
    const std::vector<int>& lib_indices,
    const std::vector<int>& pred_indices,
    const std::vector<int>& conEs,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    int totalRow,
    int totalCol,
    bool simplex,
    double theta,
    size_t threads,
    int parallel_level,
    bool cumulate,
    int style,
    int dist_metric,
    bool dist_average,
    const std::vector<int>& dir = {0}
);

/**
 * Perform Grid-based Spatially Convergent Partial Cross Mapping (SCPCM) for multiple library sizes.
 *
 * This function estimates the cross mapping and partial cross mapping between predictor variables (`xMatrix`) and response
 * variables (`yMatrix`) over a 2D spatial grid, incorporating control variables (`zMatrixs`). It supports both Simplex Projection
 * and S-Mapping, with options for parallel computation and progress tracking.
 *
 * Parameters:
 * - xMatrix: A 2D matrix of predictor variable values (spatial cross-sectional data).
 * - yMatrix: A 2D matrix of response variable values (spatial cross-sectional data).
 * - zMatrixs: A 2D matrix storing the control variables.
 * - lib_sizes: A 2D vector where the first sub-vector contains row-wise library sizes and the second sub-vector contains column-wise library sizes.
 * - lib: A vector of pairs representing the indices (row, column) of spatial units to be the library.
 * - pred: A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * - Es: A vector specifying the embedding dimensions for attractor reconstruction using `xMatrix` and control variables.
 * - taus: A vector specifying the spatial lag steps for constructing lagged state-space vectors with control variables.
 * - b: A vector specifying the numbers of nearest neighbors used for prediction.
 * - simplex: Boolean flag indicating whether to use Simplex Projection (true) or S-Mapping (false) for prediction.
 * - theta: Distance weighting parameter used for weighting neighbors in the S-Mapping prediction.
 * - threads: Number of threads to use for parallel computation.
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * - cumulate: Boolean flag indicating whether to cumulate partial correlations.
 * - style: Embedding style selector (0: includes current state, 1: excludes it).
 * - stack: Embedding arrangement selector (0: single - average lags, 1: composite - stack).  Default is 0 (average lags).
 * - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * - dist_average: Whether to average distance by the number of valid vector components.
 * - single_sig: Whether to estimate significance and confidence intervals using a single rho value.
 * - progressbar: Boolean flag indicating whether to display a progress bar during computation.
 * - dir: Direction selector for embeddings where 0 returns all directions for embeddings, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 * - win_rations: Scale the sliding window step relative to the matrix width/height to speed up state-space predictions.
 *
 * Returns:
 *    A 2D vector of results, where each row contains:
 *      - The library size.
 *      - The mean pearson cross-mapping correlation.
 *      - The statistical significance of the pearson correlation.
 *      - The lower bound of the pearson correlation confidence interval.
 *      - The upper bound of the pearson correlation confidence interval.
 *      - The mean partial cross-mapping partial correlation.
 *      - The statistical significance of the partial correlation.
 *      - The lower bound of the partial correlation confidence interval.
 *      - The upper bound of the partial correlation confidence interval.
 */
std::vector<std::vector<double>> SCPCM4Grid(
    const std::vector<std::vector<double>>& xMatrix,     // Two dimension matrix of X variable
    const std::vector<std::vector<double>>& yMatrix,     // Two dimension matrix of Y variable
    const std::vector<std::vector<double>>& zMatrixs,    // 2D matrix that stores the control variables
    const std::vector<std::vector<int>>& lib_sizes,      // Vector of library sizes to use
    const std::vector<std::pair<int, int>>& lib,         // Indices of spatial units to be the library
    const std::vector<std::pair<int, int>>& pred,        // Indices of spatial units to be predicted
    const std::vector<int>& Es,                          // Number of dimensions for the attractor reconstruction with the x and control variables
    const std::vector<int>& taus,                        // Vector specifying the spatial lag step for constructing lagged state-space vectors with control variables.
    const std::vector<int>& b,                           // Numbers of nearest neighbors to use for prediction
    bool simplex,                                        // Algorithm used for prediction; Use simplex projection if true, and s-mapping if false
    double theta,                                        // Distance weighting parameter for the local neighbours in the manifold
    int threads,                                         // Number of threads used from the global pool
    int parallel_level,                                  // Level of parallel computing: 0 for `lower`, 1 for `higher`
    bool cumulate,                                       // Whether to cumulate the partial correlations
    int style,                                           // Embedding style selector (0: includes current state, 1: excludes it)
    int stack,                                           // Embedding arrangement selector (0: single - average lags, 1: composite - stack).  Default is 0 (average lags).
    int dist_metric,                                     // Distance metric selector (1: Manhattan, 2: Euclidean)
    bool dist_average,                                   // Whether to average distance by the number of valid vector components
    bool single_sig,                                     // Whether to estimate significance and confidence intervals using a single rho value
    bool progressbar = false,                            // Whether to print the progress bar
    const std::vector<int>& dir = {0},                   // Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
    const std::vector<double>& win_ratios = {0,0}        // Scale the sliding window step relative to the matrix width/height to speed up state-space predictions.
);

/**
 * Performs Grid-based Spatially Convergent Partial Cross Mapping (SCPCM).
 *
 * Parameters:
 * - xMatrix: 2D grid of predictor variable values (row-major order)
 * - yMatrix: 2D grid of target variable values (row-major order)
 * - zMatrixs: Control variables stored as 2D grids (vector of vectors)
 * - lib_size: Number of consecutive spatial units to include in each library
 * - lib: Vector specifying library location as raster(grid) cell indices in row-major order
 * - pred: Vector specifying prediction location as raster(grid) cell indices in row-major order
 * - Es: Embedding dimensions for x and control variables
 * - taus: Spatial lag steps for x and control variables
 * - b: Numbers of nearest neighbors for prediction
 * - simplex: Use simplex projection (true) or S-mapping (false)
 * - theta: Distance weighting parameter for S-mapping
 * - threads: Number of parallel computation threads
 * - parallel_level: Level of parallel computing: 0 for `lower`, 1 for `higher`.
 * - cumulate: Enable cumulative partial correlations
 * - style: Embedding style selector (0: includes current state, 1: excludes it).
 * - stack: Embedding arrangement selector (0: single - average lags, 1: composite - stack).  Default is 0 (average lags).
 * - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean).
 * - dist_average: Whether to average distance by the number of valid vector components.
 * - single_sig: Whether to estimate significance and confidence intervals using a single rho value.
 * - progressbar: Display progress bar during computation
 * - dir: Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * Returns:
 *   2D vector containing:
 *     - Library size
 *     - Mean cross-map correlation (rho)
 *     - Rho significance
 *     - Rho lower CI
 *     - Rho upper CI
 *     - Mean partial correlation
 *     - Partial correlation significance
 *     - Partial lower CI
 *     - Partial upper CI
 */
std::vector<std::vector<double>> SCPCM4GridOneDim(
    const std::vector<std::vector<double>>& xMatrix,
    const std::vector<std::vector<double>>& yMatrix,
    const std::vector<std::vector<double>>& zMatrixs,
    const std::vector<int>& lib_sizes,
    const std::vector<int>& lib,
    const std::vector<int>& pred,
    const std::vector<int>& Es,
    const std::vector<int>& taus,
    const std::vector<int>& b,
    bool simplex,
    double theta,
    int threads,
    int parallel_level,
    bool cumulate,
    int style,
    int stack,
    int dist_metric,
    bool dist_average,
    bool single_sig,
    bool progressbar = false,
    const std::vector<int>& dir = {0}
);

#endif // SCPCM4Grid_H

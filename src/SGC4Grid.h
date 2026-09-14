#ifndef SGC4Grid_H
#define SGC4Grid_H

#include <vector>
#include <cstdint>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.

/**
 * @brief Compute directional spatial granger causality for 2D grid data using symbolic entropy measures.
 *
 * This function estimates the bidirectional spatial granger causality between two spatial variables
 * `x` and `y` observed on a 2D lattice/grid. It is based on a symbolic approximation of spatial transfer entropy,
 * which evaluates whether the spatial neighborhood of one variable improves the predictive information of the other.
 *
 * The procedure includes the following key steps:
 *
 * - **Spatial embedding construction**: Generate spatial lag vectors (`wx` and `wy`) for both variables using a neighborhood window.
 * - **Optional symbolization**: Convert raw values into discrete symbols to enhance robustness in the presence of spatial autocorrelation.
 * - **Entropy computation**: Estimate joint and marginal entropies of the original variables and their spatial embeddings.
 * - **Directional causality estimation**: Compute granger-like causal influence based on the information gain from adding the other variable's lag:
 *
 *   - Causality from X to Y:
 *     \f[
 *     SC_{x→y} = [H(y, wy) − H(wy)] − [H(y, wy, wx) − H(wy, wx)]
 *     \f]
 *
 *   - Causality from Y to X:
 *     \f[
 *     SC_{y→x} = [H(x, wx) − H(wx)] − [H(x, wx, wy) − H(wx, wy)]
 *     \f]
 *
 * - **Optional normalization**: If `normalize = true`, each causality score is scaled by its baseline information gain
 *   (i.e., the first term in each direction). This yields values bounded between -1 and 1, making them comparable across variables and scales.
 *
 * @param x         2D grid (matrix) representing variable X.
 * @param y         2D grid (matrix) representing variable Y.
 * @param dir       Direction selector for embeddings where 0 returns all directions for embeddings, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 * @param lib       A vector of pairs representing the indices (row, column) of spatial units to be the library.
 * @param pred      A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param k         Embedding neighborhood radius (e.g., k = 1 means 3×3 window).
 * @param base      Logarithm base used in entropy computation (default is 2, for bits).
 * @param symbolize Whether to discretize the data via symbolic transformation before entropy computation.
 * @param normalize Whether to normalize the causality scores to lie within [-1, 1] (default is false).
 *
 * @return A std::vector<double> of two values:
 *         - sc_x_to_y: Estimated spatial granger causality from x to y (normalized if specified).
 *         - sc_y_to_x: Estimated spatial granger causality from y to x (normalized if specified).
 */
std::vector<double> SGCSingle4Grid(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    const std::vector<int>& dir,
    const std::vector<std::pair<int, int>>& lib,
    const std::vector<std::pair<int, int>>& pred,
    size_t k,
    double base = 2,
    bool symbolize = true,
    bool normalize = false
);

/**
 * @brief Compute spatial granger causality for gridded data using spatial block bootstrap.
 *
 * This function estimates the directional spatial granger causality between two gridded variables `x` and `y`,
 * by applying a symbolic entropy-based method, and assesses the statistical significance of the causality using
 * spatial block bootstrap techniques. It calculates the causality in both directions: X → Y and Y → X.
 * Additionally, the function evaluates the significance of the estimated causality statistics by comparing them
 * to bootstrap realizations of the causality.
 *
 * The method involves the following steps:
 * - **Computation of true causality**: The function first calculates the spatial Granger causality statistic
 *   using the original data grids `x` and `y`.
 * - **Spatial block bootstrap resampling**: The grid values are resampled with spatial block bootstrapping.
 *   Each resample preserves local spatial structure and generates new bootstrap realizations of the causality statistic.
 * - **Estimation of causality for bootstrapped samples**: The causality statistic is estimated for each of the
 *   bootstrapped realizations, which involves calculating the symbolic entropy measures and their differences.
 * - **Empirical p-values**: The final p-values for both directional causality estimates (X → Y and Y → X) are
 *   derived by comparing the bootstrapped statistics with the true causality statistics.
 *
 * This approach takes into account spatial autocorrelation and allows the use of parallel processing for faster
 * bootstrap estimation. The spatial bootstrap method involves reshuffling grid cells into spatial blocks,
 * preserving local dependencies, and calculating causality for each realization.
 *
 * @param x           2D grid (matrix) of variable X.
 * @param y           2D grid (matrix) of variable Y, same size as x.
 * @param dir         Direction selector for embeddings where 0 returns all directions for embeddings, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 * @param lib         A vector of pairs representing the indices (row, column) of spatial units to be the library.
 * @param pred        A vector of pairs representing the indices (row, column) of spatial units to be predicted.
 * @param block       Vector assigning each grid cell to a spatial block for bootstrapping.
 * @param k           Neighborhood window size used for symbolization (typically 3 or 5).
 * @param threads     Number of threads to use for parallel bootstrap estimation.
 * @param boot        Number of bootstrap iterations (default: 399).
 * @param base        Base of the logarithm used in entropy computation (default: 2 for bits).
 * @param seed        Seed for the random number generator to ensure reproducibility (default: 42).
 * @param symbolize   Whether to use symbolic transformation for the grids (default: true).
 * @param normalize   Whether to normalize entropy values (optional, default: false).
 * @param progressbar Whether to show a progress bar during bootstrap computation (default: true).
 *
 * @return A vector of four values:
 *         - sc_x_to_y: Estimated spatial granger causality from x to y.
 *         - p_x_to_y: Empirical p-value for x → y based on bootstrap distribution.
 *         - sc_y_to_x: Estimated spatial granger causality from y to x.
 *         - p_y_to_x: Empirical p-value for y → x based on bootstrap distribution.
 */
std::vector<double> SGC4Grid(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& y,
    const std::vector<int>& dir,
    const std::vector<std::pair<int, int>>& lib,
    const std::vector<std::pair<int, int>>& pred,
    const std::vector<int>& block,
    int k,
    int threads,
    int boot = 399,
    double base = 2,
    unsigned long long seed = 42,
    bool symbolize = true,
    bool normalize = false,
    bool progressbar = true
);

#endif // SGC4Grid_H

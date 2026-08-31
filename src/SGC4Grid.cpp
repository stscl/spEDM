#include <vector>
#include <cstdint>
#include "CppGridUtils.h"
#include "Entropy.h"
#include "SpatialBlockBootstrap.h"
#include <RcppThread.h>

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
) {
  size_t rows = x.size();
  // size_t cols = x[0].size();

  std::vector<double> wx;
  std::vector<std::vector<double>> Ex = GenGridEmbeddings(x,1,1,1,dir);
  for (const auto& row : Ex) {
    wx.insert(wx.end(), row.begin(), row.end());
  }
  std::vector<std::vector<double>> xw = GridVec2Mat(wx,rows);

  std::vector<double> wy;
  std::vector<std::vector<double>> Ey = GenGridEmbeddings(y,1,1,1,dir);
  for (const auto& row : Ey) {
    wy.insert(wy.end(), row.begin(), row.end());
  }
  std::vector<std::vector<double>> yw = GridVec2Mat(wy,rows);

  double Hwx, Hwy, Hxwx, Hywy, Hwxwy, Hwxwyx, Hwxwyy;

  if (symbolize){
    std::vector<double> sx = GenGridSymbolization(x, lib, pred, k);
    std::vector<double> sy = GenGridSymbolization(y, lib, pred, k);
    std::vector<double> swx = GenGridSymbolization(xw, lib, pred, k);
    std::vector<double> swy = GenGridSymbolization(yw, lib, pred, k);

    std::vector<std::vector<double>> sp_series(pred.size(),std::vector<double>(4));
    for (size_t i = 0; i < pred.size(); ++i){
      sp_series[i] = {sx[i],sy[i],swx[i],swy[i]}; // 0:x 1:y 2:wx 3:wy
    }

    Hxwx = CppJoinEntropy_Disc(sp_series, {0,2}, base, true); // H(x,wx)
    Hywy = CppJoinEntropy_Disc(sp_series, {1,3}, base, true); // H(y,wy)
    Hwx = CppEntropy_Disc(swx, base, true); // H(wx)
    Hwy = CppEntropy_Disc(swy, base, true); // H(wy)
    Hwxwy = CppJoinEntropy_Disc(sp_series,{2,3}, base, true); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Disc(sp_series,{0,2,3}, base, true); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Disc(sp_series,{1,2,3}, base, true); // H(wx,wy,y)
  } else {
    std::vector<std::vector<double>> sp_series(pred.size(),std::vector<double>(4));
    std::vector<double> wx_series(pred.size());
    std::vector<double> wy_series(pred.size());
    for (size_t i = 0; i < pred.size(); ++i){
      int ri = pred[i].first;
      int ci = pred[i].second;
      sp_series[i] = {x[ri][ci],y[ri][ci],xw[ri][ci],yw[ri][ci]}; // 0:x 1:y 2:wx 3:wy
      wx_series[i] = xw[ri][ci];
      wy_series[i] = yw[ri][ci];
    }

    Hxwx = CppJoinEntropy_Cont(sp_series, {0,2}, k, base, true); // H(x,wx)
    Hywy = CppJoinEntropy_Cont(sp_series, {1,3}, k, base, true); // H(y,wy)
    Hwx = CppEntropy_Cont(wx_series, k, base, true); // H(wx)
    Hwy = CppEntropy_Cont(wy_series, k, base, true); // H(wy)
    Hwxwy = CppJoinEntropy_Cont(sp_series,{2,3}, k, base, true); // H(wx,wy)
    Hwxwyx = CppJoinEntropy_Cont(sp_series,{0,2,3}, k, base, true); // H(wx,wy,x)
    Hwxwyy = CppJoinEntropy_Cont(sp_series,{1,2,3}, k, base, true); // H(wx,wy,y)
  }

  double sc_x_to_y,sc_y_to_x;
  if (normalize){
    // transformed to fall within [-1, 1]
    sc_x_to_y = ((Hywy - Hwy) - (Hwxwyy - Hwxwy)) / (Hywy - Hwy);
    sc_y_to_x = ((Hxwx - Hwx) - (Hwxwyx - Hwxwy)) / ((Hxwx - Hwx));
  } else {
    sc_x_to_y = (Hywy - Hwy) - (Hwxwyy - Hwxwy);
    sc_y_to_x = (Hxwx - Hwx) - (Hwxwyx - Hwxwy);
  }

  return {sc_x_to_y,sc_y_to_x};
}

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
){
  // Initialize the bootstrapped realizations of the spatial granger causality statistic
  std::vector<std::vector<double>> sc_bootstraps(boot);

  const size_t rows = x.size();
  const size_t cols = x[0].size();

  // // Previous implementation may cause spurious correlations during randomization process
  // auto monte_boots = [&](int n){
  //   // Use different seed for each iteration to ensure different random samples
  //   unsigned int current_seed = seed + n;
  //   // Generate a spatial block bootstrap resample of indices
  //   std::vector<int> boot_indice = SpatialBlockBootstrap(block,current_seed);
  //   // Obtain the bootstrapped realization series
  //   std::vector<double> x_bs(boot_indice.size());
  //   std::vector<double> y_bs(boot_indice.size());
  //   for (size_t i = 0; i < boot_indice.size(); ++i){
  //     std::vector<int> cellindice = RowColFromGrid(i,cols);
  //     x_bs[i] = x[cellindice[0]][cellindice[1]];
  //     y_bs[i] = y[cellindice[0]][cellindice[1]];
  //   }
  //   std::vector<std::vector<double>> x_boot = GridVec2Mat(x_bs,static_cast<int>(rows));
  //   std::vector<std::vector<double>> y_boot = GridVec2Mat(y_bs,static_cast<int>(rows));
  //   // Estimate the bootstrapped realization of the spatial granger causality statistic
  //   sc_bootstraps[n] = SGCSingle4Grid(x_boot,y_boot,lib,pred,static_cast<size_t>(std::abs(k)),base,symbolize,normalize);
  // };

  // Prebuild 64-bit RNG pool with seed sequence
  std::vector<std::mt19937_64> rng_pool(boot);
  for (int i = 0; i < boot; ++i) {
    std::seed_seq seq{static_cast<uint64_t>(seed), static_cast<uint64_t>(i)};
    rng_pool[i] = std::mt19937_64(seq);
  }

  auto monte_boots = [&](int n){
    // Use prebuilt rng instance
    std::vector<int> boot_indice = SpatialBlockBootstrapRNG(block, rng_pool[n]);
    // Obtain the bootstrapped realization series
    std::vector<double> x_bs(boot_indice.size());
    std::vector<double> y_bs(boot_indice.size());
    for (size_t i = 0; i < boot_indice.size(); ++i){
      std::vector<int> cellindice = RowColFromGrid(i,cols);
      x_bs[i] = x[cellindice[0]][cellindice[1]];
      y_bs[i] = y[cellindice[0]][cellindice[1]];
    }
    std::vector<std::vector<double>> x_boot = GridVec2Mat(x_bs,static_cast<int>(rows));
    std::vector<std::vector<double>> y_boot = GridVec2Mat(y_bs,static_cast<int>(rows));
    // Estimate the bootstrapped realization of the spatial granger causality statistic
    sc_bootstraps[n] = SGCSingle4Grid(x_boot,y_boot,dir,lib,pred,static_cast<size_t>(std::abs(k)),base,symbolize,normalize);
  };

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Parallel computation with or without a progress bar
  if (progressbar) {
    RcppThread::ProgressBar bar(boot, 1);
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
      bar++;
    }, threads_sizet);
  } else {
    RcppThread::parallelFor(0, boot, [&](int i) {
      monte_boots(i);
    }, threads_sizet);
  }

  // The "true" spatial granger causality statistic
  std::vector<double> sc = SGCSingle4Grid(x,y,dir,lib,pred,static_cast<size_t>(std::abs(k)),base,symbolize,normalize);
  double scx = sc[0];
  double scy = sc[1];
  // Compute the estimated bootstrap p–value
  double b_xy = 0;
  double b_yx = 0;
  for (size_t i = 0; i < sc_bootstraps.size(); ++i){
    if (sc_bootstraps[i][0] > scx){
      b_xy += 1;
    }
    if (sc_bootstraps[i][1] > scy) {
      b_yx += 1;
    }
  }

  return {scx,b_xy / boot,scy,b_yx / boot};
}

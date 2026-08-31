#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <tuple>
#include <unordered_map>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionalCardinality.h"
#include "PatternCausality.h"
#include <RcppThread.h>

/*
 * Evaluates prediction performance of different combinations of embedding dimensions, number of nearest neighbors and
 * tau values for grid data using simplex projection forecasting.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - E: A vector of embedding dimensions to evaluate.
 *   - b: A vector of nearest neighbors to use for prediction.
 *   - tau: A vector of spatial lag steps for constructing lagged state-space vectors.
 *   - style: Embedding style selector (0: includes current state, 1: excludes it).  Default is 1 (excludes current state).
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
 *   - dir: Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * Returns:
 *   A 2D vector where each row contains [E, b, tau, rho, mae, rmse] for a given embedding dimension.
 */
std::vector<std::vector<double>> Simplex4Grid(const std::vector<std::vector<double>>& source,
                                              const std::vector<std::vector<double>>& target,
                                              const std::vector<int>& lib_indices,
                                              const std::vector<int>& pred_indices,
                                              const std::vector<int>& E,
                                              const std::vector<int>& b,
                                              const std::vector<int>& tau,
                                              int style = 1,
                                              int dist_metric = 2,
                                              bool dist_average = true,
                                              int threads = 8,
                                              const std::vector<int>& dir = {0}) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const int numRows = target.size();
  const int numCols = target[0].size();

  // Flatten target matrix
  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols);
  for (const auto& row : target) {
    vec_std.insert(vec_std.end(), row.begin(), row.end());
  }

  // Unique sorted embedding dimensions, neighbor values, and tau values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, b, tau) combinations
  std::vector<std::tuple<int, int, int>> unique_EbTau;
  for (int e : Es) {
    for (int bb : bs) {
      if (bb < e) continue;

      for (int t : taus)
        unique_EbTau.emplace_back(e, bb, t);
    }
  }

  std::vector<std::vector<double>> result(unique_EbTau.size(), std::vector<double>(6));

  // Parallel loop over combinations
  RcppThread::parallelFor(0, unique_EbTau.size(), [&](size_t i) {
    const int cur_E = std::get<0>(unique_EbTau[i]);
    const int cur_b = std::get<1>(unique_EbTau[i]);
    const int cur_tau = std::get<2>(unique_EbTau[i]);
    // auto [cur_E, cur_b, cur_tau] = unique_EbTau[i]; // C++17 structured binding

    // Generate embedding
    std::vector<std::vector<double>> embeddings = GenGridEmbeddings(source, cur_E, cur_tau, style, dir);

    // Evaluate performance
    std::vector<double> metrics = SimplexBehavior(embeddings, vec_std, lib_indices, pred_indices, cur_b, dist_metric, dist_average);

    // Store results
    result[i][0] = cur_E;
    result[i][1] = cur_b;
    result[i][2] = cur_tau;
    result[i][3] = metrics[0];
    result[i][4] = metrics[1];
    result[i][5] = metrics[2];
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different combinations of embedding dimensions, number of nearest neighbors and
 * tau values for grid data using simplex projection forecasting (composite embeddings version).
 */
std::vector<std::vector<double>> Simplex4GridCom(const std::vector<std::vector<double>>& source,
                                                 const std::vector<std::vector<double>>& target,
                                                 const std::vector<int>& lib_indices,
                                                 const std::vector<int>& pred_indices,
                                                 const std::vector<int>& E,
                                                 const std::vector<int>& b,
                                                 const std::vector<int>& tau,
                                                 int style = 1,
                                                 int dist_metric = 2,
                                                 bool dist_average = true,
                                                 int threads = 8,
                                                 const std::vector<int>& dir = {0}) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const int numRows = target.size();
  const int numCols = target[0].size();

  // Flatten target matrix
  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols);
  for (const auto& row : target) {
    vec_std.insert(vec_std.end(), row.begin(), row.end());
  }

  // Unique sorted embedding dimensions, neighbor values, and tau values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, b, tau) combinations
  std::vector<std::tuple<int, int, int>> unique_EbTau;
  for (int e : Es) {
    for (int bb : bs) {
      if (bb < e) continue;

      for (int t : taus)
        unique_EbTau.emplace_back(e, bb, t);
    }
  }

  std::vector<std::vector<double>> result(unique_EbTau.size(), std::vector<double>(6));

  // Parallel loop over combinations
  RcppThread::parallelFor(0, unique_EbTau.size(), [&](size_t i) {
    const int cur_E = std::get<0>(unique_EbTau[i]);
    const int cur_b = std::get<1>(unique_EbTau[i]);
    const int cur_tau = std::get<2>(unique_EbTau[i]);
    // auto [cur_E, cur_b, cur_tau] = unique_EbTau[i]; // C++17 structured binding

    // Generate embedding
    std::vector<std::vector<std::vector<double>>> embeddings = GenGridEmbeddingsCom(source, cur_E, cur_tau, style, dir);

    // Evaluate performance
    std::vector<double> metrics = SimplexBehavior(embeddings, vec_std, lib_indices, pred_indices, cur_b, dist_metric, dist_average);

    // Store results
    result[i][0] = cur_E;
    result[i][1] = cur_b;
    result[i][2] = cur_tau;
    result[i][3] = metrics[0];
    result[i][4] = metrics[1];
    result[i][5] = metrics[2];
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for grid data
 * using the S-mapping method.
 *
 * Parameters:
 *   - source: A matrix to be embedded.
 *   - target: A matrix to be predicted.
 *   - lib_indices: A vector of indices indicating the library (training) set.
 *   - pred_indices: A vector of indices indicating the prediction set.
 *   - theta: A vector of weighting parameters for distance calculation in SMap.
 *   - E: The embedding dimension to evaluate. Default is 3.
 *   - tau: The spatial lag step for constructing lagged state-space vectors. Default is 1.
 *   - b: Number of nearest neighbors to use for prediction. Default is 4.
 *   - style: Embedding style selector (0: includes current state, 1: excludes it).  Default is 1 (excludes current state).
 *   - dist_metric: Distance metric selector (1: Manhattan, 2: Euclidean). Default is 2 (Euclidean).
 *   - dist_average: Whether to average distance by the number of valid vector components. Default is true.
 *   - threads: Number of threads used from the global pool. Default is 8.
 *   - dir: Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * Returns:
 *   A 2D vector where each row contains [theta, rho, mae, rmse] for a given theta value.
 */
std::vector<std::vector<double>> SMap4Grid(const std::vector<std::vector<double>>& source,
                                           const std::vector<std::vector<double>>& target,
                                           const std::vector<int>& lib_indices,
                                           const std::vector<int>& pred_indices,
                                           const std::vector<double>& theta,
                                           int E = 3,
                                           int tau = 1,
                                           int b = 4,
                                           int style = 1,
                                           int dist_metric = 2,
                                           bool dist_average = true,
                                           int threads = 8,
                                           const std::vector<int>& dir = {0}) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const int numRows = target.size();
  const int numCols = target[0].size();

  // Flatten target matrix
  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols);
  for (const auto& row : target) {
    vec_std.insert(vec_std.end(), row.begin(), row.end());
  }

  // Generate embedding once
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(source, E, tau, style, dir);

  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    std::vector<double> metrics = SMapBehavior(embeddings, vec_std, lib_indices, pred_indices, b, theta[i], dist_average, dist_metric);

    result[i][0] = theta[i];
    result[i][1] = metrics[0];
    result[i][2] = metrics[1];
    result[i][3] = metrics[2];
  }, threads_sizet);

  return result;
}

/*
 * Evaluates prediction performance of different theta parameters for grid data
 * using the S-mapping method (composite embeddings version).
 */
std::vector<std::vector<double>> SMap4GridCom(const std::vector<std::vector<double>>& source,
                                              const std::vector<std::vector<double>>& target,
                                              const std::vector<int>& lib_indices,
                                              const std::vector<int>& pred_indices,
                                              const std::vector<double>& theta,
                                              int E = 3,
                                              int tau = 1,
                                              int b = 4,
                                              int style = 1,
                                              int dist_metric = 2,
                                              bool dist_average = true,
                                              int threads = 8,
                                              const std::vector<int>& dir = {0}) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  const int numRows = target.size();
  const int numCols = target[0].size();

  // Flatten target matrix
  std::vector<double> vec_std;
  vec_std.reserve(numRows * numCols);
  for (const auto& row : target) {
    vec_std.insert(vec_std.end(), row.begin(), row.end());
  }

  // Generate embedding once
  std::vector<std::vector<std::vector<double>>> embeddings = GenGridEmbeddingsCom(source, E, tau, style, dir);

  std::vector<std::vector<double>> result(theta.size(), std::vector<double>(4));

  RcppThread::parallelFor(0, theta.size(), [&](size_t i) {
    std::vector<double> metrics = SMapBehavior(embeddings, vec_std, lib_indices, pred_indices, b, theta[i], dist_average, dist_metric);

    result[i][0] = theta[i];
    result[i][1] = metrics[0];
    result[i][2] = metrics[1];
    result[i][3] = metrics[2];
  }, threads_sizet);

  return result;
}

/**
 * @brief Evaluate intersectional cardinality (IC) for spatial grid data.
 *
 * This function computes the intersectional cardinality between the k-nearest neighbors
 * of grid-embedded source and target spatial variables, across a range of embedding dimensions (E),
 * neighborhood sizes (b) and spatial lag step (tau). The result is an AUC (Area Under the Curve)
 * score for each (E, tau) pair that quantifies the directional similarity or interaction between
 * the spatial fields.
 *
 * The method constructs delay-like embeddings over grid cells using spatial neighborhoods,
 * filters out invalid prediction locations (e.g., with all NaN values), computes nearest neighbors
 * in embedding space, and calculates the cardinality of overlapping neighbors. These overlaps are
 * then evaluated using a CMC-based statistical test (via AUC).
 *
 * Supports both single-threaded and parallel execution using `RcppThread`.
 *
 * @param source 2D spatial variable (grid) used as the source for embedding.
 * @param target 2D spatial variable (grid) used as the target for embedding.
 * @param lib_indices Indices of spatial locations used as the library set (training).
 * @param pred_indices Indices of spatial locations used as the prediction set (evaluation).
 * @param E Vector of spatial embedding dimensions to evaluate (e.g., neighborhood sizes).
 * @param b Vector of neighbor counts (k) used to compute IC.
 * @param tau Vector of spatial embedding spacing (lag). Determines step between spatial lag.
 * @param exclude Number of nearest neighbors to exclude in IC computation.
 * @param style Embedding style selector (0: includes current state, 1: excludes it).
 * @param dist_metric Distance metric selector (1: Manhattan, 2: Euclidean).
 * @param threads Maximum number of threads to use.
 * @param parallel_level If > 0, enables parallel evaluation of b for each E.
 * @param dir Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * @return A matrix of size (|E| × |b| × |tau|) × 5 with rows: [E, b, tau, AUC, P-value]
 */
std::vector<std::vector<double>> IC4Grid(const std::vector<std::vector<double>>& source,
                                         const std::vector<std::vector<double>>& target,
                                         const std::vector<size_t>& lib_indices,
                                         const std::vector<size_t>& pred_indices,
                                         const std::vector<int>& E,
                                         const std::vector<int>& b,
                                         const std::vector<int>& tau,
                                         int exclude = 0,
                                         int style = 1,
                                         int dist_metric = 2,
                                         int threads = 8,
                                         int parallel_level = 0,
                                         const std::vector<int>& dir = {0}) {
  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Remove duplicates from E and b
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, tau) combinations
  std::vector<std::pair<int, int>> unique_ETau;
  unique_ETau.reserve(Es.size() * taus.size());
  for (int e : Es)
    for (int t : taus)
      unique_ETau.emplace_back(e, t);

  std::vector<std::vector<double>> result(unique_ETau.size() * bs.size(),
                                          std::vector<double>(5));

  size_t max_num_neighbors = 0;
  if (!bs.empty()) {
    max_num_neighbors = static_cast<size_t>(bs.back() + exclude);
  }

  if (parallel_level == 0){
    for (size_t i = 0; i < unique_ETau.size(); ++i) {
      const int Ei = unique_ETau[i].first;
      const int taui = unique_ETau[i].second;

      // Generate embeddings
      auto embedding_x = GenGridEmbeddings(source, Ei, taui, style, dir);
      auto embedding_y = GenGridEmbeddings(target, Ei, taui, style, dir);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<size_t> valid_pred;
      for (size_t idx : pred_indices) {
        if (idx < 0 || idx >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
      bool L1norm = (dist_metric == 1);

      // // Precompute neighbors (The earlier implementation based on a serial version)
      // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true),lib_indices,max_num_neighbors);
      // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true),lib_indices,max_num_neighbors);

      // Precompute neighbors (parallel computation)
      auto nx = CppMatKNNeighbors(embedding_x, lib_indices, max_num_neighbors, threads_sizet, L1norm);
      auto ny = CppMatKNNeighbors(embedding_y, lib_indices, max_num_neighbors, threads_sizet, L1norm);

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      for (size_t j = 0; j < bs.size(); ++j){
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionalCardinalitySingle(
          nx,ny,lib_indices.size(),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,0
        );

        std::vector<double> cs = {0,1};
        if (!res.empty())  cs = CppCMCTest(res[0].Intersection,">");

        result[j + bs.size() * i][0] = Ei;     // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = taui;   // tau
        result[j + bs.size() * i][3] = cs[0];  // AUC
        result[j + bs.size() * i][4] = cs[1];  // P value
      }
    }
  } else {
    for (size_t i = 0; i < unique_ETau.size(); ++i) {
      const int Ei = unique_ETau[i].first;
      const int taui = unique_ETau[i].second;

      // Generate embeddings
      auto embedding_x = GenGridEmbeddings(source, Ei, taui, style, dir);
      auto embedding_y = GenGridEmbeddings(target, Ei, taui, style, dir);

      // Filter valid prediction points (exclude those with all NaN values)
      std::vector<size_t> valid_pred;
      for (size_t idx : pred_indices) {
        if (idx < 0 || idx >= embedding_x.size()) continue;

        bool x_nan = std::all_of(embedding_x[idx].begin(), embedding_x[idx].end(),
                                 [](double v) { return std::isnan(v); });
        bool y_nan = std::all_of(embedding_y[idx].begin(), embedding_y[idx].end(),
                                 [](double v) { return std::isnan(v); });
        if (!x_nan && !y_nan) valid_pred.push_back(idx);
      }

      // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
      bool L1norm = (dist_metric == 1);

      // // Precompute neighbors (The earlier implementation based on a serial version)
      // auto nx = CppDistSortedIndice(CppMatDistance(embedding_x, L1norm, true),lib_indices,max_num_neighbors);
      // auto ny = CppDistSortedIndice(CppMatDistance(embedding_y, L1norm, true),lib_indices,max_num_neighbors);

      // Precompute neighbors (parallel computation)
      auto nx = CppMatKNNeighbors(embedding_x, lib_indices, max_num_neighbors, threads_sizet, L1norm);
      auto ny = CppMatKNNeighbors(embedding_y, lib_indices, max_num_neighbors, threads_sizet, L1norm);

      // Parameter initialization
      const size_t n_excluded_sizet = static_cast<size_t>(exclude);

      RcppThread::parallelFor(0, bs.size(), [&](size_t j) {
        const size_t k = static_cast<size_t>(bs[j]);

        // run cross mapping
        std::vector<IntersectionRes> res = IntersectionalCardinalitySingle(
          nx,ny,lib_indices.size(),lib_indices,valid_pred,k,n_excluded_sizet,threads_sizet,1
        );

        std::vector<double> cs = {0,1};
        if (!res.empty())  cs = CppCMCTest(res[0].Intersection,">");

        result[j + bs.size() * i][0] = Ei;     // E
        result[j + bs.size() * i][1] = bs[j];  // k
        result[j + bs.size() * i][2] = taui;   // tau
        result[j + bs.size() * i][3] = cs[0];  // AUC
        result[j + bs.size() * i][4] = cs[1];  // P value
      }, threads_sizet);
    }
  }

  return result;
}

/**
 * @brief Search for optimal embedding parameters in geographical pattern causality analysis on regular grids.
 *
 * @param source
 *   A 2D matrix (vector of vectors) representing the source spatiotemporal field
 *   on a regular grid. Each inner vector corresponds to a spatial location’s time series.
 *
 * @param target
 *   A 2D matrix representing the target field, aligned with `source` in both
 *   spatial index and temporal length.
 *
 * @param lib_indices
 *   Indices used as the "library set" for neighbor searching and pattern estimation.
 *   Typically corresponds to time indices used for model fitting.
 *
 * @param pred_indices
 *   Indices used as the "prediction set", where causality is evaluated.
 *
 * @param E
 *   Candidate embedding dimensions. Each E defines the number of delay coordinates
 *   used to reconstruct state vectors in the grid embedding.
 *
 * @param b
 *   Candidate neighbor counts. Each b specifies how many nearest neighbors are used
 *   in the pattern similarity / causality computation.
 *
 * @param tau
 *   Candidate spatial lags (lag steps) used when constructing delay embeddings.
 *
 * @param style
 *   Embedding style for grid reconstruction:
 *   (See GenGridEmbeddings() for supported styles.)
 *
 * @param zero_tolerance
 *   Threshold for treating small distances or weights as zero.
 *   Useful for avoiding numerical instability.
 *
 * @param dist_metric
 *   Distance metric used for neighbor searching:
 *     - 1: Manhattan distance
 *     - 2: Euclidean distance (default)
 *
 * @param relative
 *   Whether to normalize distances relative to local scale (relative = true),
 *   or treat them as absolute values.
 *
 * @param weighted
 *   Whether to weight causal strength.
 *
 * @param threads
 *   Maximum number of threads to use for parallel computation.
 *   If negative, absolute value is used. If larger than hardware limit,
 *   it is automatically capped.
 *
 * @param parallel_level
 *   Controls the parallel level of computation.
 *
 * @param dir
 *   Direction selector for embeddings where 0 returns all directions, 1–8 correspond to NW, N, NE, W, E, SW, S, SE, and multiple directions can be combined (e.g., {1,2,3} for NW, N, NE).
 *
 * @return
 *   A 2D matrix where each row corresponds to one (E, b, tau) parameter triplet:
 *
 *     [0] E
 *     [1] b
 *     [2] tau
 *     [3] TotalPos   — Strength of positive-pattern causality
 *     [4] TotalNeg   — Strength of negative-pattern causality
 *     [5] TotalDark  — Strength of ambiguous / non-directional causality
 *
 *   This table enables searching for optimal embedding parameters in
 *   geographical Pattern Causality studies.
 */
std::vector<std::vector<double>> PC4Grid(const std::vector<std::vector<double>>& source,
                                         const std::vector<std::vector<double>>& target,
                                         const std::vector<size_t>& lib_indices,
                                         const std::vector<size_t>& pred_indices,
                                         const std::vector<int>& E,
                                         const std::vector<int>& b,
                                         const std::vector<int>& tau,
                                         int style = 1,
                                         int zero_tolerance = 0,
                                         int dist_metric = 2,
                                         bool relative = true,
                                         bool weighted = true,
                                         int threads = 8,
                                         int parallel_level = 0,
                                         const std::vector<int>& dir = {0}) {
  // Unique sorted embedding dimensions, neighbor values, and tau values
  std::vector<int> Es = E;
  std::sort(Es.begin(), Es.end());
  Es.erase(std::unique(Es.begin(), Es.end()), Es.end());

  std::vector<int> bs = b;
  std::sort(bs.begin(), bs.end());
  bs.erase(std::unique(bs.begin(), bs.end()), bs.end());

  std::vector<int> taus = tau;
  std::sort(taus.begin(), taus.end());
  taus.erase(std::unique(taus.begin(), taus.end()), taus.end());

  // Generate unique (E, b, tau) combinations
  std::vector<std::tuple<int, int, int>> unique_EbTau;
  for (int e : Es) {
    for (int bb : bs) {
      if (bb < e) continue;

      for (int t : taus)
        unique_EbTau.emplace_back(e, bb, t);
    }
  }

  // Prepare for data slicing
  std::vector<size_t> selected_indices;
  selected_indices.reserve(lib_indices.size() + pred_indices.size());
  for (size_t i = 0; i < lib_indices.size(); ++i)
    selected_indices.push_back(lib_indices[i]);
  for (size_t i = 0; i < pred_indices.size(); ++i)
    selected_indices.push_back(pred_indices[i]);
  std::sort(selected_indices.begin(), selected_indices.end());
  selected_indices.erase(
    std::unique(selected_indices.begin(), selected_indices.end()),
    selected_indices.end()
  );

  // Check if full set is used
  bool use_subset = (selected_indices.size() < target.size() * target[0].size());
  
  std::vector<size_t> lib_sub, pred_sub;
  lib_sub.reserve(lib_indices.size());
  pred_sub.reserve(pred_indices.size());

  if (use_subset) {
    std::unordered_map<size_t, size_t> index_map;
    index_map.reserve(selected_indices.size());

    for (size_t i = 0; i < selected_indices.size(); ++i) {
      index_map[selected_indices[i]] = i;
    }

    // --- Remap lib indices ---
    for (size_t i = 0; i < lib_indices.size(); ++i) {
      lib_sub.push_back(index_map[lib_indices[i]]);
    }

    // --- Remap pred indices ---
    for (size_t i = 0; i < pred_indices.size(); ++i) {
      pred_sub.push_back(index_map[pred_indices[i]]);
    }
  }

  std::vector<std::vector<double>> result(unique_EbTau.size(), std::vector<double>(6));

  if (parallel_level == 0) {
    for (size_t i = 0; i < unique_EbTau.size(); ++i) {
      const int Ei   = std::get<0>(unique_EbTau[i]);
      const int bi   = std::get<1>(unique_EbTau[i]);
      const int taui = std::get<2>(unique_EbTau[i]);
      // auto [Ei, bi, taui] = unique_EbTau[i]; // C++17 structured binding

      auto Mx = GenGridEmbeddings(source, Ei, taui, style, dir);
      auto My = GenGridEmbeddings(target, Ei, taui, style, dir);
      
      PatternCausalityRes res;

      if (!use_subset) {
        // --- Full data: no slicing needed ---
        res = PatternCausality(
          Mx, My, lib_indices, pred_indices, bi, zero_tolerance,
          dist_metric, relative, weighted, threads);
      } else {
        // --- Slice Mx and My ---
        std::vector<std::vector<double>> Mx_sub;
        std::vector<std::vector<double>> My_sub;

        Mx_sub.reserve(selected_indices.size());
        My_sub.reserve(selected_indices.size());

        for (size_t i = 0; i < selected_indices.size(); ++i) {
          size_t idx = selected_indices[i];
          Mx_sub.push_back(Mx[idx]);
          My_sub.push_back(My[idx]);
        }

        // --- Compute pattern causality on subset ---
        res = PatternCausality(
          Mx_sub, My_sub, lib_sub, pred_sub, bi, zero_tolerance,
          dist_metric, relative, weighted, threads);
      }

      result[i][0] = Ei;
      result[i][1] = bi;
      result[i][2] = taui;
      result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
      result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
      result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
    }
  } else {
    // Configure threads
    size_t threads_sizet = static_cast<size_t>(std::abs(threads));
    threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

    RcppThread::parallelFor(0, unique_EbTau.size(), [&](size_t i) {
      const int Ei   = std::get<0>(unique_EbTau[i]);
      const int bi   = std::get<1>(unique_EbTau[i]);
      const int taui = std::get<2>(unique_EbTau[i]);
      // auto [Ei, bi, taui] = unique_EbTau[i]; // C++17 structured binding

      auto Mx = GenGridEmbeddings(source, Ei, taui, style, dir);
      auto My = GenGridEmbeddings(target, Ei, taui, style, dir);

      PatternCausalityRes res;

      if (!use_subset) {
        // --- Full data: no slicing needed ---
        res = PatternCausality(
          Mx, My, lib_indices, pred_indices, bi, zero_tolerance,
          dist_metric, relative, weighted, 1);
      } else {
        // --- Slice Mx and My ---
        std::vector<std::vector<double>> Mx_sub;
        std::vector<std::vector<double>> My_sub;

        Mx_sub.reserve(selected_indices.size());
        My_sub.reserve(selected_indices.size());

        for (size_t i = 0; i < selected_indices.size(); ++i) {
          size_t idx = selected_indices[i];
          Mx_sub.push_back(Mx[idx]);
          My_sub.push_back(My[idx]);
        }

        // --- Compute pattern causality on subset ---
        res = PatternCausality(
          Mx_sub, My_sub, lib_sub, pred_sub, bi, zero_tolerance,
          dist_metric, relative, weighted, 1);
      }

      result[i][0] = Ei;
      result[i][1] = bi;
      result[i][2] = taui;
      result[i][3] = std::isnan(res.TotalPos) ? 0.0 : res.TotalPos;
      result[i][4] = std::isnan(res.TotalNeg) ? 0.0 : res.TotalNeg;
      result[i][5] = std::isnan(res.TotalDark) ? 0.0 : res.TotalDark;
    }, threads_sizet);
  }

  return result;
}

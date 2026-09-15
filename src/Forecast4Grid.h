#ifndef Forecast4Grid_H
#define Forecast4Grid_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <tuple>
#include "CppGridUtils.h"
#include "SimplexProjection.h"
#include "SMap.h"
#include "IntersectionalCardinality.h"
#include "PatternCausality.h"
// Note: <RcppThread.h> is intentionally excluded from this header to avoid
//       unnecessary Rcpp dependencies and potential header inclusion order
//       issues (e.g., R.h being included before Rcpp headers). It should only
//       be included in the corresponding .cpp implementation file.

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
                                              const std::vector<int>& dir = {0});

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
                                                 const std::vector<int>& dir = {0});

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
                                           const std::vector<int>& dir = {0});

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
                                              const std::vector<int>& dir = {0});

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
                                         const std::vector<int>& dir = {0});

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
                                         const std::vector<int>& dir = {0});

#endif // Forecast4Grid_H

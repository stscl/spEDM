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
#include <RcppThread.h>

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
) {
  // Extract row-wise and column-wise library sizes
  const int lib_size_row = lib_sizes[0];
  const int lib_size_col = lib_sizes[1];

  // Determine the marked libsize
  const int libsize = row_size_mark ? lib_size_row : lib_size_col;

  // Computate step of sliding windows
  const int r_step = std::max(1, static_cast<int>(std::round(1.0 + win_ratios[0] * lib_size_row)));
  const int c_step = std::max(1, static_cast<int>(std::round(1.0 + win_ratios[1] * lib_size_col)));

  // Precompute valid (r, c) pairs
  std::vector<std::pair<int, int>> valid_indices;
  // // Slide the window by 1 row and 1 column
  // for (int r = 1; r <= totalRow - lib_size_row + 1; ++r) {
  //   for (int c = 1; c <= totalCol - lib_size_col + 1; ++c) {
  //     valid_indices.emplace_back(r, c);
  //   }
  // }
  // Slide the window by flexible row and column
  for (int r = 1; r <= totalRow - lib_size_row + 1; r += r_step) {
    for (int c = 1; c <= totalCol - lib_size_col + 1; c += c_step) {
      valid_indices.emplace_back(r, c);
    }
  }

  // // Initialize the result container with the same size as valid_indices
  // std::vector<std::pair<int, double>> x_xmap_y;
  // x_xmap_y.resize(valid_indices.size());

  // Preallocate the result vector to avoid out-of-bounds access
  std::vector<std::pair<int, double>> x_xmap_y(valid_indices.size());

  // Unified processing logic
  auto process = [&](size_t idx) {
    const int r = valid_indices[idx].first;
    const int c = valid_indices[idx].second;

    // Initialize library indices and count the number of the nan value together
    std::vector<int> lib_indices;

    for (int ii = r; ii < r + lib_size_row; ++ii) {
      for (int jj = c; jj < c + lib_size_col; ++jj) {
        int index = (ii - 1) * totalCol + (jj - 1);
        if (possible_lib_indices[index]) {
          lib_indices.emplace_back(index);
        }
      }
    }

    double rho = std::numeric_limits<double>::quiet_NaN();
    // Run cross map and store results
    if (simplex) {
      rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b, dist_metric, dist_average);
    } else {
      rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta, dist_metric, dist_average);
    }
    x_xmap_y[idx] = std::make_pair(libsize, rho);
  };

  // Parallel coordination
  if (parallel_level == 0) {
    RcppThread::parallelFor(0, valid_indices.size(), process, threads);
  } else {
    for (size_t i = 0; i < valid_indices.size(); ++i) {
      process(i);
    }
  }

  return x_xmap_y;
}

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
) {
  // Extract row-wise and column-wise library sizes
  const int lib_size_row = lib_sizes[0];
  const int lib_size_col = lib_sizes[1];

  // Determine the marked libsize
  const int libsize = row_size_mark ? lib_size_row : lib_size_col;

  // Computate step of sliding windows
  const int r_step = std::max(1, static_cast<int>(std::round(1.0 + win_ratios[0] * lib_size_row)));
  const int c_step = std::max(1, static_cast<int>(std::round(1.0 + win_ratios[1] * lib_size_col)));

  // Precompute valid (r, c) pairs
  std::vector<std::pair<int, int>> valid_indices;
  // // Slide the window by 1 row and 1 column
  // for (int r = 1; r <= totalRow - lib_size_row + 1; ++r) {
  //   for (int c = 1; c <= totalCol - lib_size_col + 1; ++c) {
  //     valid_indices.emplace_back(r, c);
  //   }
  // }
  // Slide the window by flexible row and column
  for (int r = 1; r <= totalRow - lib_size_row + 1; r += r_step) {
    for (int c = 1; c <= totalCol - lib_size_col + 1; c += c_step) {
      valid_indices.emplace_back(r, c);
    }
  }

  // // Initialize the result container with the same size as valid_indices
  // std::vector<std::pair<int, double>> x_xmap_y;
  // x_xmap_y.resize(valid_indices.size());

  // Preallocate the result vector to avoid out-of-bounds access
  std::vector<std::pair<int, double>> x_xmap_y(valid_indices.size());

  // Unified processing logic
  auto process = [&](size_t idx) {
    const int r = valid_indices[idx].first;
    const int c = valid_indices[idx].second;

    // Initialize library indices and count the number of the nan value together
    std::vector<int> lib_indices;

    for (int ii = r; ii < r + lib_size_row; ++ii) {
      for (int jj = c; jj < c + lib_size_col; ++jj) {
        int index = (ii - 1) * totalCol + (jj - 1);
        if (possible_lib_indices[index]) {
          lib_indices.emplace_back(index);
        }
      }
    }

    double rho = std::numeric_limits<double>::quiet_NaN();
    // Run cross map and store results
    if (simplex) {
      rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b, dist_metric, dist_average);
    } else {
      rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta, dist_metric, dist_average);
    }
    x_xmap_y[idx] = std::make_pair(libsize, rho);
  };

  // Parallel coordination
  if (parallel_level == 0) {
    RcppThread::parallelFor(0, valid_indices.size(), process, threads);
  } else {
    for (size_t i = 0; i < valid_indices.size(); ++i) {
      process(i);
    }
  }

  return x_xmap_y;
}

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
) {
  int max_lib_size = lib_indices.size();

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<std::pair<int, double>> x_xmap_y;

    double rho = std::numeric_limits<double>::quiet_NaN();
    // Run cross map and store results
    if (simplex) {
      rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b, dist_metric, dist_average);
    } else {
      rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta, dist_metric, dist_average);
    }
    x_xmap_y.emplace_back(lib_size, rho);
    return x_xmap_y;
  } else if (parallel_level == 0){
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<std::pair<int, double>> x_xmap_y(valid_lib_indices.size());

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
      double rho = std::numeric_limits<double>::quiet_NaN();
      // Run cross map and store results
      if (simplex) {
        rho = SimplexProjection(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, dist_metric, dist_average);
      } else {
        rho = SMap(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, theta, dist_metric, dist_average);
      }

      std::pair<int, double> result(lib_size, rho); // Store the product of row and column library sizes
      x_xmap_y[i] = result;
    }, threads);

    return x_xmap_y;
  } else {
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<std::pair<int, double>> x_xmap_y(valid_lib_indices.size());

    // Iterate through precomputed valid_lib_indices
    for (size_t i = 0; i < valid_lib_indices.size(); ++i){
      double rho = std::numeric_limits<double>::quiet_NaN();
      // Run cross map and store results
      if (simplex) {
        rho = SimplexProjection(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, dist_metric, dist_average);
      } else {
        rho = SMap(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, theta, dist_metric, dist_average);
      }

      std::pair<int, double> result(lib_size, rho); // Store the product of row and column library sizes
      x_xmap_y[i] = result;
    }

    return x_xmap_y;
  }
}

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
) {
  int max_lib_size = lib_indices.size();

  // No possible library variation if using all vectors
  if (lib_size == max_lib_size) {
    std::vector<std::pair<int, double>> x_xmap_y;

    double rho = std::numeric_limits<double>::quiet_NaN();
    // Run cross map and store results
    if (simplex) {
      rho = SimplexProjection(xEmbedings, yPred, lib_indices, pred_indices, b, dist_metric, dist_average);
    } else {
      rho = SMap(xEmbedings, yPred, lib_indices, pred_indices, b, theta, dist_metric, dist_average);
    }
    x_xmap_y.emplace_back(lib_size, rho);
    return x_xmap_y;
  } else if (parallel_level == 0){
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<std::pair<int, double>> x_xmap_y(valid_lib_indices.size());

    // Perform the operations using RcppThread
    RcppThread::parallelFor(0, valid_lib_indices.size(), [&](size_t i) {
      double rho = std::numeric_limits<double>::quiet_NaN();
      // Run cross map and store results
      if (simplex) {
        rho = SimplexProjection(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, dist_metric, dist_average);
      } else {
        rho = SMap(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, theta, dist_metric, dist_average);
      }

      std::pair<int, double> result(lib_size, rho); // Store the product of row and column library sizes
      x_xmap_y[i] = result;
    }, threads);

    return x_xmap_y;
  } else {
    // Precompute valid indices for the library
    std::vector<std::vector<int>> valid_lib_indices;
    for (int start_lib = 0; start_lib < max_lib_size; ++start_lib) {
      std::vector<int> local_lib_indices;
      // Loop around to beginning of lib indices
      if (start_lib + lib_size > max_lib_size) {
        for (int i = start_lib; i < max_lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
        int num_vectors_remaining = lib_size - (max_lib_size - start_lib);
        for (int i = 0; i < num_vectors_remaining; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      } else {
        for (int i = start_lib; i < start_lib + lib_size; ++i) {
          local_lib_indices.emplace_back(lib_indices[i]);
        }
      }
      valid_lib_indices.emplace_back(local_lib_indices);
    }

    // Preallocate the result vector to avoid out-of-bounds access
    std::vector<std::pair<int, double>> x_xmap_y(valid_lib_indices.size());

    // Iterate through precomputed valid_lib_indices
    for (size_t i = 0; i < valid_lib_indices.size(); ++i){
      double rho = std::numeric_limits<double>::quiet_NaN();
      // Run cross map and store results
      if (simplex) {
        rho = SimplexProjection(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, dist_metric, dist_average);
      } else {
        rho = SMap(xEmbedings, yPred, valid_lib_indices[i], pred_indices, b, theta, dist_metric, dist_average);
      }

      std::pair<int, double> result(lib_size, rho); // Store the product of row and column library sizes
      x_xmap_y[i] = result;
    }

    return x_xmap_y;
  }
}

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
) {
  // If b is not provided correctly, default it to E + 2
  if (b <= 0) {
    b = E + 2;
  }

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> Emb2D;
  std::vector<std::vector<std::vector<double>>> Emb3D;
  if (stack == 0){
    Emb2D = GenGridEmbeddings(xMatrix, E, tau, style, dir);
  } else {
    Emb3D = GenGridEmbeddingsCom(xMatrix, E, tau, style, dir);
  }

  // Ensure the maximum value does not exceed totalRow or totalCol
  int max_lib_size_row = totalRow;
  int max_lib_size_col = totalCol;

  // Extract row-wise and column-wise library sizes
  std::vector<int> row_lib_sizes = lib_sizes[0];
  std::vector<int> col_lib_sizes = lib_sizes[1];

  // Transform to ensure no size exceeds max_lib_size_row or max_lib_size_col
  std::transform(row_lib_sizes.begin(), row_lib_sizes.end(), row_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size_row); });
  std::transform(col_lib_sizes.begin(), col_lib_sizes.end(), col_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size_col); });

  // Remove duplicates in row-wise and column-wise library sizes
  std::sort(row_lib_sizes.begin(), row_lib_sizes.end());
  row_lib_sizes.erase(std::unique(row_lib_sizes.begin(), row_lib_sizes.end()), row_lib_sizes.end());
  std::sort(col_lib_sizes.begin(), col_lib_sizes.end());
  col_lib_sizes.erase(std::unique(col_lib_sizes.begin(), col_lib_sizes.end()), col_lib_sizes.end());

  // // Generate unique pairs of row-wise and column-wise library sizes
  // std::vector<std::pair<int, int>> unique_lib_size_pairs;
  // for (int row_size : row_lib_sizes) {
  //   for (int col_size : col_lib_sizes) {
  //     unique_lib_size_pairs.emplace_back(row_size, col_size);
  //   }
  // }

  // Generate unique pairs of row-wise and column-wise library sizes
  std::vector<std::pair<int, int>> unique_lib_size_pairs;

  // Determine which library size vector is longer
  int row_size_count = row_lib_sizes.size();
  int col_size_count = col_lib_sizes.size();
  int min_size = std::min(row_size_count, col_size_count);
  // int max_size = std::max(row_size_count, col_size_count);

  // Fill unique_lib_size_pairs based on the shorter vector
  for (int i = 0; i < min_size; ++i) {
    unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes[i]);
  }

  bool row_size_mark = true;
  // Handle the excess elements for the longer vector
  if (row_size_count > col_size_count) {
    for (int i = min_size; i < row_size_count; ++i) {
      unique_lib_size_pairs.emplace_back(row_lib_sizes[i], col_lib_sizes.back()); // Pair with the max value of col_lib_sizes
    }
  }

  if (row_size_count < col_size_count) {
    for (int i = min_size; i < col_size_count; ++i) {
      unique_lib_size_pairs.emplace_back(row_lib_sizes.back(), col_lib_sizes[i]); // Pair with the max value of row_lib_sizes
    }
    row_size_mark = false;
  }

  // Set library indices
  std::vector<bool> lib_indices(totalRow * totalCol, false);
  for (const auto& l : lib) {
    lib_indices[LocateGridIndices(l.first + 1, l.second + 1, totalRow, totalCol)] = true;
  }

  // Set predict indices
  std::vector<int> pred_indices;
  for (const auto& p : pred) {
    pred_indices.push_back(LocateGridIndices(p.first + 1, p.second + 1, totalRow, totalCol));
  }

  // Local results for each library
  std::vector<std::vector<std::pair<int, double>>> local_results(unique_lib_size_pairs.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
      for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        if (stack == 0){
          local_results[i] = GCCMSingle4Grid(
            Emb2D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        } else {
          local_results[i] = GCCMSingle4Grid(
            Emb3D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        }
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_size_pairs.size(); ++i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        if (stack == 0){
          local_results[i] = GCCMSingle4Grid(
            Emb2D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        } else {
          local_results[i] = GCCMSingle4Grid(
            Emb3D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        }
      }
    }
  } else {
    // Perform the operations using RcppThread
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_size_pairs.size(), 1);
      RcppThread::parallelFor(0, unique_lib_size_pairs.size(), [&](size_t i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        if (stack == 0){
          local_results[i] = GCCMSingle4Grid(
            Emb2D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        } else {
          local_results[i] = GCCMSingle4Grid(
            Emb3D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        }
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_size_pairs.size(), [&](size_t i) {
        int lib_size_row = unique_lib_size_pairs[i].first;
        int lib_size_col = unique_lib_size_pairs[i].second;
        if (stack == 0){
          local_results[i] = GCCMSingle4Grid(
            Emb2D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        } else {
          local_results[i] = GCCMSingle4Grid(
            Emb3D,
            yPred,
            {lib_size_row, lib_size_col},
            lib_indices,
            pred_indices,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            row_size_mark,
            dist_metric,
            dist_average,
            win_ratios);
        }
      }, threads_sizet);
    }
  }

  // Initialize the result container
  std::vector<std::pair<int, double>> x_xmap_y;

  // Merge all local results into the final result
  for (const auto& local_result : local_results) {
    x_xmap_y.insert(x_xmap_y.end(), local_result.begin(), local_result.end());
  }

  // Group by the first int (library size) and compute the mean
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  size_t n = pred.size();
  std::vector<std::vector<double>> final_results;

  if (single_sig) {
    // Calculate significance and confidence intervals using the mean of rho vector only.
    for (const auto& group : grouped_results) {
      double mean_value = CppMean(group.second, true);
      final_results.push_back({static_cast<double>(group.first), mean_value});
    }

    // Calculate significance and confidence interval for each result
    for (size_t i = 0; i < final_results.size(); ++i) {
      double rho = final_results[i][1];
      double significance = CppCorSignificance(rho, n);
      std::vector<double> confidence_interval = CppCorConfidence(rho, n);

      final_results[i].push_back(significance);
      final_results[i].push_back(confidence_interval[0]);
      final_results[i].push_back(confidence_interval[1]);
    }
  } else {
    // Compute significance and confidence intervals directly from grouped correlation vectors
    for (const auto& group : grouped_results) {
      // Calculate the mean correlation coefficient from the group
      double mean_value = CppMean(group.second, true);

      // Compute significance (p-value) using the vector of correlations directly
      double significance = CppMeanCorSignificance(group.second, n);

      // Compute confidence interval using the vector of correlations directly
      std::vector<double> confidence_interval = CppMeanCorConfidence(group.second, n);

      // Store results: group ID, mean correlation, p-value, lower CI, upper CI
      final_results.push_back({
        static_cast<double>(group.first),
        mean_value,
        significance,
        confidence_interval[0],
        confidence_interval[1]
      });
    }
  }

  return final_results;
}

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
) {
  // If b is not provided correctly, default it to E + 2
  if (b <= 0) {
    b = E + 2;
  }

  // Configure threads
  size_t threads_sizet = static_cast<size_t>(std::abs(threads));
  threads_sizet = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), threads_sizet);

  // Get the dimensions of the xMatrix
  int totalRow = xMatrix.size();
  int totalCol = xMatrix[0].size();

  // Flatten yMatrix into a 1D array (row-major order)
  std::vector<double> yPred;
  for (const auto& row : yMatrix) {
    yPred.insert(yPred.end(), row.begin(), row.end());
  }

  // Generate embeddings for xMatrix
  std::vector<std::vector<double>> Emb2D;
  std::vector<std::vector<std::vector<double>>> Emb3D;
  if (stack == 0){
    Emb2D = GenGridEmbeddings(xMatrix, E, tau, style, dir);
  } else {
    Emb3D = GenGridEmbeddingsCom(xMatrix, E, tau, style, dir);
  }

  int max_lib_size = static_cast<int>(lib.size()); // Maximum lib size

  std::vector<int> unique_lib_sizes(lib_sizes.begin(), lib_sizes.end());

  // Transform to ensure no size exceeds max_lib_size
  std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
                 [&](int size) { return std::min(size, max_lib_size); });

  // Ensure the minimum value in unique_lib_sizes is b (uncomment this section if required)
  // std::transform(unique_lib_sizes.begin(), unique_lib_sizes.end(), unique_lib_sizes.begin(),
  //                [&](int size) { return std::max(size, b); });

  // Remove duplicates
  std::sort(unique_lib_sizes.begin(), unique_lib_sizes.end());
  unique_lib_sizes.erase(std::unique(unique_lib_sizes.begin(), unique_lib_sizes.end()), unique_lib_sizes.end());

  // Local results for each library
  std::vector<std::vector<std::pair<int, double>>> local_results(unique_lib_sizes.size());

  if (parallel_level == 0){
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        if (stack == 0){
          local_results[i] = GCCMSingle4GridOneDim(
            Emb2D,
            yPred,
            unique_lib_sizes[i],
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        } else {
          local_results[i] = GCCMSingle4GridOneDim(
            Emb3D,
            yPred,
            unique_lib_sizes[i],
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        }
        bar++;
      }
    } else {
      for (size_t i = 0; i < unique_lib_sizes.size(); ++i) {
        if (stack == 0){
          local_results[i] = GCCMSingle4GridOneDim(
            Emb2D,
            yPred,
            unique_lib_sizes[i],
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        } else {
          local_results[i] = GCCMSingle4GridOneDim(
            Emb3D,
            yPred,
            unique_lib_sizes[i],
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        }
      }
    }
  } else {
    // Iterate over each library size
    if (progressbar) {
      RcppThread::ProgressBar bar(unique_lib_sizes.size(), 1);
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        if (stack == 0){
          local_results[i] = GCCMSingle4GridOneDim(
            Emb2D,
            yPred,
            lib_size,
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        } else {
          local_results[i] = GCCMSingle4GridOneDim(
            Emb3D,
            yPred,
            lib_size,
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        }
        bar++;
      }, threads_sizet);
    } else {
      RcppThread::parallelFor(0, unique_lib_sizes.size(), [&](size_t i) {
        int lib_size = unique_lib_sizes[i];
        if (stack == 0){
          local_results[i] = GCCMSingle4GridOneDim(
            Emb2D,
            yPred,
            lib_size,
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        } else {
          local_results[i] = GCCMSingle4GridOneDim(
            Emb3D,
            yPred,
            lib_size,
            lib,
            pred,
            totalRow,
            totalCol,
            b,
            simplex,
            theta,
            threads_sizet,
            parallel_level,
            dist_metric,
            dist_average
          );
        }
      }, threads_sizet);
    }
  }

  // Initialize the result container
  std::vector<std::pair<int, double>> x_xmap_y;

  // Merge all local results into the final result
  for (const auto& local_result : local_results) {
    x_xmap_y.insert(x_xmap_y.end(), local_result.begin(), local_result.end());
  }

  // Group by the first int (library size) and compute the mean
  std::map<int, std::vector<double>> grouped_results;
  for (const auto& result : x_xmap_y) {
    grouped_results[result.first].push_back(result.second);
  }

  size_t n = pred.size();
  std::vector<std::vector<double>> final_results;

  if (single_sig) {
    // Calculate significance and confidence intervals using the mean of rho vector only.
    for (const auto& group : grouped_results) {
      double mean_value = CppMean(group.second, true);
      final_results.push_back({static_cast<double>(group.first), mean_value});
    }

    // Calculate significance and confidence interval for each result
    for (size_t i = 0; i < final_results.size(); ++i) {
      double rho = final_results[i][1];
      double significance = CppCorSignificance(rho, n);
      std::vector<double> confidence_interval = CppCorConfidence(rho, n);

      final_results[i].push_back(significance);
      final_results[i].push_back(confidence_interval[0]);
      final_results[i].push_back(confidence_interval[1]);
    }
  } else {
    // Compute significance and confidence intervals directly from grouped correlation vectors
    for (const auto& group : grouped_results) {
      // Calculate the mean correlation coefficient from the group
      double mean_value = CppMean(group.second, true);

      // Compute significance (p-value) using the vector of correlations directly
      double significance = CppMeanCorSignificance(group.second, n);

      // Compute confidence interval using the vector of correlations directly
      std::vector<double> confidence_interval = CppMeanCorConfidence(group.second, n);

      // Store results: group ID, mean correlation, p-value, lower CI, upper CI
      final_results.push_back({
        static_cast<double>(group.first),
        mean_value,
        significance,
        confidence_interval[0],
        confidence_interval[1]
      });
    }
  }

  return final_results;
}

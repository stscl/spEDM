#include <vector>
#include <cmath>
#include <string>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include "CppStats.h"
#include "CppGridUtils.h"
#include "Forecast4Grid.h"
#include "MultiViewEmbedding.h"
#include "GCCM4Grid.h"
#include "SCPCM4Grid.h"
#include "CrossMappingCardinality.h"
#include "PatternCausality.h"
#include "FalseNearestNeighbors.h"
#include "SLM4Grid.h"
#include "SGC4Grid.h"
// 'Rcpp.h' should not be included and correct to include only 'RcppArmadillo.h'.
// #include <Rcpp.h>
#include <RcppArmadillo.h>

// Wrapper function to convert row and column indices to 1-based linear indices
// [[Rcpp::export(rng = false)]]
int RcppLocateGridIndices(int curRow, int curCol,
                          int totalRow, int totalCol){
  int indices = LocateGridIndices(curRow,curCol,totalRow,totalCol);
  return indices + 1;
};

// Wrapper function to convert a row-major linear index to row and column indices
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppRowColFromGrid(int cellNum, int totalCol){
  std::vector<int> result = RowColFromGrid(cellNum - 1, totalCol);
  for (int& val : result) {
    val += 1;
  }
  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(result);
};

// Wrapper function to calculate lagged values for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppLaggedVal4Grid(const Rcpp::NumericMatrix& mat, int lagNum = 1,
                                       const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0–8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // Call the CppLaggedVal4Grid function
  std::vector<std::vector<double>> laggedMat = CppLaggedVal4Grid(cppMat, lagNum, dir_std);

  // Convert the result back to Rcpp::NumericMatrix
  int laggedRows = laggedMat.size();
  int laggedCols = laggedMat[0].size();
  Rcpp::NumericMatrix result(laggedRows, laggedCols);

  for (int r = 0; r < laggedRows; ++r) {
    for (int c = 0; c < laggedCols; ++c) {
      result(r, c) = laggedMat[r][c];
    }
  }

  return result;
}

// Wrapper function to generate embeddings for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGenGridEmbeddings(const Rcpp::NumericMatrix& mat,
                                          int E = 3, int tau = 1, int style = 1,
                                          const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // Call the GenGridEmbeddings function
  std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, E, tau, style, dir_std);

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  int rows = embeddings.size();
  int cols = embeddings[0].size();
  Rcpp::NumericMatrix result(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      result(i, j) = embeddings[i][j];
    }
  }

  return result;
}

// Wrapper function to generate composite grid embeddings.
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGenGridEmbeddingsCom(const Rcpp::NumericMatrix& mat,
                                    int E = 3,
                                    int tau = 1,
                                    int style = 1,
                                    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert R matrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // Call the core C++ embedding function
  std::vector<std::vector<std::vector<double>>> embeddings =
    GenGridEmbeddingsCom(cppMat, E, tau, style, dir_std);

  // Convert the 3D std::vector into an Rcpp::List of NumericMatrix
  Rcpp::List result;

  for (const auto& layer : embeddings) {
    if (layer.empty()) continue; // Skip empty layers (if all NaN subsets were removed)
    int rows = layer.size();
    int cols = layer[0].size();
    Rcpp::NumericMatrix mat_layer(rows, cols);

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        mat_layer(i, j) = layer[i][j];
      }
    }

    result.push_back(mat_layer);
  }

  return result;
}

// Wrapper function to generate neighbors for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGenGridNeighbors(const Rcpp::NumericMatrix& mat,
                                const Rcpp::IntegerMatrix& lib,
                                int k = 8) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<int> lib_std(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std[i] = lib(i, 0) - 1;
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std[i] = LocateGridIndices(lib(i, 0), lib(i, 1), numRows, numCols);
    }
  }

  // Generate neighbors
  std::vector<std::vector<int>> neighbors = GenGridNeighbors(
    cppMat, lib_std, static_cast<size_t>(std::abs(k))
  );

  // Convert neighbors to Rcpp::List with 1-based indexing
  int n = neighbors.size();
  Rcpp::List result(n);
  for (int i = 0; i < n; ++i) {
    std::vector<int> neighbor_i = neighbors[i];
    for (auto& idx : neighbor_i) {
      idx += 1; // convert to 1-based index
    }
    result[i] = Rcpp::IntegerVector(neighbor_i.begin(), neighbor_i.end());
  }

  return result;
}

// Wrapper function to implement a symbolic transformation of a univariate spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppGenGridSymbolization(const Rcpp::NumericMatrix& mat,
                                             const Rcpp::IntegerMatrix& lib,
                                             const Rcpp::IntegerMatrix& pred,
                                             int k = 8) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<std::pair<int, int>> lib_std;
  lib_std.reserve(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
      lib_std.emplace_back(rowcolnum[0], rowcolnum[1]);
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      lib_std.emplace_back(lib(i, 0) - 1, lib(i, 1) - 1);
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<std::pair<int, int>> pred_std;
  pred_std.reserve(pred.nrow());

  if (pred_dim == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
      pred_std.emplace_back(rowcolnum[0], rowcolnum[1]);
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      pred_std.emplace_back(pred(i, 0) - 1, pred(i, 1) - 1);
    }
  }

  // Call the GenGridSymbolization function
  std::vector<double> symbolmap = GenGridSymbolization(
    cppMat, lib_std, pred_std, static_cast<size_t>(std::abs(k))
  );

  // Convert the result back to Rcpp::NumericVector
  return Rcpp::wrap(symbolmap);
}

// Wrapper function to partition spatial units in spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::IntegerVector RcppDivideGrid(const Rcpp::NumericMatrix& mat,
                                   int b = 9, int shape = 3) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the CppDivideGrid function
  std::vector<int> blocks = CppDivideGrid(cppMat, b, shape);

  // Convert the result back to Rcpp::IntegerVector
  return Rcpp::wrap(blocks);
}

// Wrapper function to perform univariate Spatial Logistic Map for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSLMUni4Grid(
    const Rcpp::NumericMatrix& mat,
    int k = 4,
    int step = 20,
    double alpha = 0.77,
    double noise_level = 0.0,
    double escape_threshold = 1e10,
    unsigned long long random_seed = 42
) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Call the core function
  std::vector<std::vector<double>> result = SLMUni4Grid(cppMat, k, step, alpha, noise_level, 
                                                        escape_threshold, random_seed);

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result.size());
  int n_cols = result[0].size();
  Rcpp::NumericMatrix out(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out(i, j) = result[i][j];
    }
  }

  return out;
}

// Wrapper function to perform bivariate Spatial Logistic Map for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSLMBi4Grid(
    const Rcpp::NumericMatrix& mat1,
    const Rcpp::NumericMatrix& mat2,
    int k = 4,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double beta_xy = 0.05,
    double beta_yx = 0.4,
    int interact = 0,
    double noise_level = 0.0,
    double escape_threshold = 1e10,
    unsigned long long random_seed = 42
) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat1.nrow();
  int numCols = mat1.ncol();
  std::vector<std::vector<double>> cppMat1(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> cppMat2(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat1[r][c] = mat1(r, c);
      cppMat2[r][c] = mat2(r, c);
    }
  }

  // Call the core function
  std::vector<std::vector<std::vector<double>>> result = SLMBi4Grid(
    cppMat1, cppMat2, k, step, alpha_x, alpha_y, beta_xy, beta_yx,
    interact, noise_level, escape_threshold, random_seed
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result[0].size());
  int n_cols = result[0][1].size();
  Rcpp::NumericMatrix out_x(n_rows, n_cols);
  Rcpp::NumericMatrix out_y(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out_x(i, j) = result[0][i][j];
      out_y(i, j) = result[1][i][j];
    }
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y
  );

  return out;
}

// Wrapper function to perform trivariate Spatial Logistic Map for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppSLMTri4Grid(
    const Rcpp::NumericMatrix& mat1,
    const Rcpp::NumericMatrix& mat2,
    const Rcpp::NumericMatrix& mat3,
    int k = 4,
    int step = 20,
    double alpha_x = 0.625,
    double alpha_y = 0.77,
    double alpha_z = 0.55,
    double beta_xy = 0.05,
    double beta_xz = 0.05,
    double beta_yx = 0.4,
    double beta_yz = 0.4,
    double beta_zx = 0.65,
    double beta_zy = 0.65,
    int interact = 0,
    double noise_level = 0.0,
    double escape_threshold = 1e10,
    unsigned long long random_seed = 42
) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat1.nrow();
  int numCols = mat1.ncol();
  std::vector<std::vector<double>> cppMat1(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> cppMat2(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> cppMat3(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat1[r][c] = mat1(r, c);
      cppMat2[r][c] = mat2(r, c);
      cppMat3[r][c] = mat3(r, c);
    }
  }

  // Call the core function
  std::vector<std::vector<std::vector<double>>> result = SLMTri4Grid(
    cppMat1, cppMat2, cppMat3,
    k, step, alpha_x, alpha_y, alpha_z,
    beta_xy, beta_xz, beta_yx, beta_yz, beta_zx, beta_zy,
    interact, noise_level, escape_threshold, random_seed
  );

  // Create NumericMatrix with rows = number of spatial units, cols = number of steps+1
  int n_rows = static_cast<int>(result[0].size());
  int n_cols = result[0][1].size();
  Rcpp::NumericMatrix out_x(n_rows, n_cols);
  Rcpp::NumericMatrix out_y(n_rows, n_cols);
  Rcpp::NumericMatrix out_z(n_rows, n_cols);

  // Copy data into NumericMatrix
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      out_x(i, j) = result[0][i][j];
      out_y(i, j) = result[1][i][j];
      out_z(i, j) = result[2][i][j];
    }
  }

  // Wrap results into an Rcpp::List
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("x") = out_x,
    Rcpp::Named("y") = out_y,
    Rcpp::Named("z") = out_z
  );

  return out;
}

// Wrapper function to perform FNN for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppFNN4Grid(
    const Rcpp::NumericMatrix& mat,
    const Rcpp::NumericVector& rt,
    const Rcpp::NumericVector& eps,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    int tau = 1,
    int style = 1,
    int stack = 0,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0,
    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)){
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = mat.nrow();
  int numCols = mat.ncol();
  std::vector<std::vector<double>> cppMat(numRows, std::vector<double>(numCols));
  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      cppMat[r][c] = mat(r, c);
    }
  }

  // Convert Rcpp NumericVector to std::vector<double>
  std::vector<double> rt_std = Rcpp::as<std::vector<double>>(rt);
  std::vector<double> eps_std = Rcpp::as<std::vector<double>>(eps);

  // Convert Rcpp IntegerMatrix to std::vector<int>
  int n_libcol = lib.ncol();
  int n_predcol = pred.ncol();

  std::vector<size_t> lib_std;
  lib_std.reserve(lib.nrow());
  if (n_libcol == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(cppMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_std.push_back(static_cast<size_t>(lib(i,0) - 1));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int rowLibIndice = lib(i,0);
      int colLibIndice = lib(i,1);
      if (!std::isnan(cppMat[rowLibIndice-1][colLibIndice-1])){
        lib_std.push_back(static_cast<size_t>(LocateGridIndices(rowLibIndice, colLibIndice, numRows, numCols)));
      }
    }
  }

  std::vector<size_t> pred_std;
  pred_std.reserve(pred.nrow());
  if (n_predcol == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(cppMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_std.push_back(static_cast<size_t>(pred(i,0) - 1));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int rowPredIndice = pred(i,0);
      int colPredIndice = pred(i,1);
      if (!std::isnan(cppMat[rowPredIndice-1][colPredIndice-1])){
        pred_std.push_back(static_cast<size_t>(LocateGridIndices(rowPredIndice, colPredIndice, numRows, numCols)));
      }
    }
  }

  // Use L1 norm (Manhattan distance) if dist_metric == 1, else use L2 norm
  bool L1norm = (dist_metric == 1);

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0–8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // Generate embeddings and perform FNN for spatial grid data
  std::vector<double> E_std = Rcpp::as<std::vector<double>>(E);
  int max_E = CppMax(E_std, true);

  std::vector<double> fnn;
  if (stack == 0){
    std::vector<std::vector<double>> embeddings = GenGridEmbeddings(cppMat, max_E, tau, style, dir_std);
    fnn = CppFNN(embeddings,lib_std,pred_std,rt_std,eps_std,L1norm,threads,parallel_level);
  } else {
    std::vector<std::vector<std::vector<double>>> embeddings = GenGridEmbeddingsCom(cppMat, max_E, tau, style, dir_std);
    fnn = CppFNN(embeddings,lib_std,pred_std,rt_std,eps_std,L1norm,threads,parallel_level);
  }

  // Convert the result back to Rcpp::NumericVector and set names as "E:1", "E:2", ..., "E:n"
  Rcpp::NumericVector result = Rcpp::wrap(fnn);
  Rcpp::CharacterVector resnames(result.size());
  for (int i = 0; i < result.size(); ++i) {
    resnames[i] = "E:" + std::to_string(i + 1);
  }
  result.names() = resnames;

  // Terminal-friendly hint (one-time, non-intrusive)
  Rcpp::Rcout << "[fnn] Input E values exceeding max embeddable dimension were truncated, and values < 2 were clamped to 2.\n"
              << "[fnn] Max embedding dimension E_max is auto-computed, with results returned for dimensions 1 through E_max.\n"
              << "[fnn] Output 'E:i' (where i = 1 to E_max-1) corresponds to the comparison between dimension i and i+1.\n";

  return result;
}

// Wrapper function to perform parameter selection of simplex forecasting for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSimplex4Grid(const Rcpp::NumericMatrix& source,
                                     const Rcpp::NumericMatrix& target,
                                     const Rcpp::IntegerMatrix& lib,
                                     const Rcpp::IntegerMatrix& pred,
                                     const Rcpp::IntegerVector& E,
                                     const Rcpp::IntegerVector& b,
                                     const Rcpp::IntegerVector& tau,
                                     int style = 1,
                                     int stack = 0,
                                     int dist_metric = 2,
                                     bool dist_average = true,
                                     int threads = 8,
                                     const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  lib_indices.reserve(lib.nrow());
  std::vector<int> pred_indices;
  pred_indices.reserve(pred.nrow());

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(lib(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        lib_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(pred(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        pred_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> b_std = Rcpp::as<std::vector<int>>(b);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  std::vector<std::vector<double>> res_std;
  if (stack == 0){
    res_std = Simplex4Grid(
      sourceMat,
      targetMat,
      lib_indices,
      pred_indices,
      E_std,
      b_std,
      tau_std,
      style,
      dist_metric,
      dist_average,
      threads,
      dir_std);
  } else {
    res_std = Simplex4GridCom(
      sourceMat,
      targetMat,
      lib_indices,
      pred_indices,
      E_std,
      b_std,
      tau_std,
      style,
      dist_metric,
      dist_average,
      threads,
      dir_std);
  }

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "tau", "rho", "mae", "rmse");
  return result;
}

// Wrapper function to perform parameter selection of s-mapping for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSMap4Grid(const Rcpp::NumericMatrix& source,
                                  const Rcpp::NumericMatrix& target,
                                  const Rcpp::IntegerMatrix& lib,
                                  const Rcpp::IntegerMatrix& pred,
                                  const Rcpp::NumericVector& theta,
                                  int E = 3,
                                  int tau = 1,
                                  int b = 5,
                                  int style = 1,
                                  int stack = 0,
                                  int dist_metric = 2,
                                  bool dist_average = true,
                                  int threads = 8,
                                  const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  lib_indices.reserve(lib.nrow());
  std::vector<int> pred_indices;
  pred_indices.reserve(pred.nrow());

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(lib(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        lib_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(pred(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        pred_indices.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
    }
  }

  // Convert Rcpp::NumericMatrix to std::vector<double>
  std::vector<double> theta_std = Rcpp::as<std::vector<double>>(theta);

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  std::vector<std::vector<double>> res_std;
  if (stack == 0){
    res_std = SMap4Grid(
      sourceMat,
      targetMat,
      lib_indices,
      pred_indices,
      theta_std,
      E,
      tau,
      b,
      style,
      dist_metric,
      dist_average,
      threads,
      dir_std);
  } else {
    res_std = SMap4GridCom(
      sourceMat,
      targetMat,
      lib_indices,
      pred_indices,
      theta_std,
      E,
      tau,
      b,
      style,
      dist_metric,
      dist_average,
      threads,
      dir_std);
  }

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("theta", "rho", "mae", "rmse");
  return result;
}

// Wrapper function to perform multiview forecasting for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppMultiView4Grid(const Rcpp::NumericMatrix& xMatrix,
                                       const Rcpp::NumericMatrix& yMatrix,
                                       const Rcpp::IntegerMatrix& lib,
                                       const Rcpp::IntegerMatrix& pred,
                                       int E = 3,
                                       int tau = 1,
                                       int b = 5,
                                       int top = 5,
                                       int nvar = 3,
                                       int style = 1,
                                       int stack = 0,
                                       int dist_metric = 2,
                                       int dist_average = true,
                                       int threads = 8,
                                       const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)){
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(numRows, std::vector<double>(numCols));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert yMatrix to std::vector<double>
  std::vector<double> target(numRows * numCols);

  if (numCols == 1){
    for (int i = 0; i < numRows; ++i) {
      target[i] =  yMatrix(i,0);
    }
  } else {
    target = GridMat2Vec(yMatrix_cpp);
  }

  // Initialize lib_indices and pred_indices
  std::vector<int> lib_indices;
  lib_indices.reserve(lib.nrow());
  std::vector<int> pred_indices;
  pred_indices.reserve(pred.nrow());

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(yMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(lib(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      int currow = lib(i,0);
      int curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      int cellindice = LocateGridIndices(currow, curcol, numRows, numCols);
      if (!std::isnan(target[cellindice])){
        lib_indices.push_back(cellindice);
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(pred(i,0)-1);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      int currow = pred(i,0);
      int curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      int cellindice = LocateGridIndices(currow, curcol, numRows, numCols);
      if (!std::isnan(target[cellindice])){
        pred_indices.push_back(cellindice);
      }
    }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  int num_row = xMatrix.nrow();
  int num_var = xMatrix.ncol();

  //  if top <= 0, we choose to apply the heuristic of k (sqrt(m))
  int k;
  if (top <= 0){
    double m = CppCombine(num_var*E,nvar) - CppCombine(num_var*(E - 1),nvar);
    k = std::floor(std::sqrt(m));
  } else {
    k = top;
  }

  std::vector<double> res;
  if (stack == 0) { // ---- CASE 1: standard 2D embedding (stack == 0) ----
    std::vector<std::vector<double>> vec_std;
    vec_std.reserve(num_row);  // preallocate number of rows
    // Initialize rows with empty vectors
    for (int i = 0; i < num_row; ++i)
      vec_std.emplace_back();  // create num_row empty rows
    // vec_std.resize(num_row); or using resize

    for (int n = 0; n < num_var; ++n) {
      // Extract nth column
      std::vector<double> univec(num_row);

      // Copy the nth column from the matrix to the vector
      for (int i = 0; i < num_row; ++i) {
        univec[i] = xMatrix(i, n);  // Access element at (i, n)
      }
      std::vector<std::vector<double>> unimat = GridVec2Mat(univec,numRows);

      // Generate the embedding:
      std::vector<std::vector<double>> vectors = GenGridEmbeddings(unimat,E,tau,style,dir_std);

      // Append columns from embedding into existing rows (column-wise stacking)
      for (int row = 0; row < num_row; ++row) {
        vec_std[row].insert(vec_std[row].end(),
                            std::make_move_iterator(vectors[row].begin()),
                            std::make_move_iterator(vectors[row].end()));
      }
    }

    // // Filter invalid (NaN) columns
    // std::vector<size_t> validColumns;
    // for (size_t col = 0; col < vec_std[0].size(); ++col) {
    //   bool isAllNaN = true;
    //   for (size_t row = 0; row < vec_std.size(); ++row) {
    //     if (!std::isnan(vec_std[row][col])) { isAllNaN = false; break; }
    //   }
    //   if (!isAllNaN) validColumns.push_back(col);
    // }
    //
    // if (validColumns.size() != vec_std[0].size()) {
    //   std::vector<std::vector<double>> filteredEmbeddings;
    //   filteredEmbeddings.reserve(vec_std.size());
    //   for (const auto& row : vec_std) {
    //     std::vector<double> filteredRow;
    //     for (size_t col : validColumns) filteredRow.push_back(row[col]);
    //     filteredEmbeddings.push_back(std::move(filteredRow));
    //   }
    //   vec_std.swap(filteredEmbeddings);
    // }

    // Perform multi-view embedding (2D version)
    res = MultiViewEmbedding(
      vec_std, target, lib_indices, pred_indices,
      b, k, dist_metric, dist_average, threads);

  } else {  // ---- CASE 2: stacked 3D embedding (stack != 0) ----
    // stacked_vec will store embeddings stacked along the column dimension.
    // Structure: stacked_vec[embedding_dimension][row][col]
    std::vector<std::vector<std::vector<double>>> stacked_vec;

    // Initialize stacked_vec based on the first variable's shape
    {
      std::vector<double> univec(num_row);
      for (int i = 0; i < num_row; ++i) univec[i] = xMatrix(i, 0);
      std::vector<std::vector<double>> unimat = GridVec2Mat(univec,numRows);
      auto embedding = GenGridEmbeddingsCom(unimat, E, tau, style, dir_std);

      // Initialize stacked_vec with correct shape
      stacked_vec = std::move(embedding);
    }

    // Start from variable index 1 since index 0 is initialized
    for (int n = 1; n < num_var; ++n) {
      // Extract variable column
      std::vector<double> univec(num_row);
      for (int i = 0; i < num_row; ++i) univec[i] = xMatrix(i, n);
      std::vector<std::vector<double>> unimat = GridVec2Mat(univec,numRows);

      // Get embedding for this variable
      auto embedding = GenGridEmbeddingsCom(unimat, E, tau, style, dir_std);

      // Append each embedding block column-wise
      for (size_t j = 0; j < stacked_vec.size(); ++j) {
        // Each stacked_vec[j] and embedding[j] must share same row count
        for (size_t r = 0; r < stacked_vec[j].size(); ++r) {
          // Append columns from embedding[j][r] to stacked_vec[j][r]
          stacked_vec[j][r].insert(
              stacked_vec[j][r].end(),
              std::make_move_iterator(embedding[j][r].begin()),
              std::make_move_iterator(embedding[j][r].end())
          );
        }
      }
    }

    // Perform multi-view embedding (3D version)
    res = MultiViewEmbedding(
      stacked_vec, target, lib_indices, pred_indices,
      b, k, dist_metric, dist_average, threads);
  }

  // Initialize a NumericMatrix with the given dimensions
  Rcpp::NumericMatrix resmat(numRows, numCols);

  // Fill the matrix with elements from res (assuming row-major order)
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      resmat(i, j) = res[i * numCols + j];  // Access element in row-major order
    }
  }

  return resmat;
}

// Wrapper function to perform parameter selection of intersectional cardinality for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppIC4Grid(const Rcpp::NumericMatrix& source,
                                const Rcpp::NumericMatrix& target,
                                const Rcpp::IntegerMatrix& lib,
                                const Rcpp::IntegerMatrix& pred,
                                const Rcpp::IntegerVector& E,
                                const Rcpp::IntegerVector& b,
                                const Rcpp::IntegerVector& tau,
                                int exclude = 0,
                                int style = 1,
                                int dist_metric = 2,
                                int threads = 8,
                                int parallel_level = 0,
                                const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<size_t> lib_indices;
  lib_indices.reserve(lib.nrow());
  std::vector<size_t> pred_indices;
  pred_indices.reserve(pred.nrow());

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
        lib_indices.push_back(static_cast<size_t>(lib(i,0)-1));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        lib_indices.push_back(static_cast<size_t>(LocateGridIndices(currow, curcol, numRows, numCols)));
      }
    }
  }

  if (pred_col == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_indices.push_back(static_cast<size_t>(pred(i,0)-1));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
        pred_indices.push_back(static_cast<size_t>(LocateGridIndices(currow, curcol, numRows, numCols)));
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  // Check the validity of the neignbor numbers
  std::vector<int> b_std;
  for (int i = 0; i < b.size(); ++i){
    if (b[i] > static_cast<int>(lib_indices.size())) {
      Rcpp::stop("Neighbor numbers count out of acceptable range at position %d (value: %d)", i + 1, b[i]);
    }
    b_std.push_back(b[i]);
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  std::vector<std::vector<double>> res_std = IC4Grid(
    sourceMat,
    targetMat,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau_std,
    exclude,
    style,
    dist_metric,
    threads,
    parallel_level,
    dir_std);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "tau", "CausalScore", "Significance");
  return result;
}

// Wrapper function to compute pattern causality for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppPC4Grid(const Rcpp::NumericMatrix& source,
                                const Rcpp::NumericMatrix& target,
                                const Rcpp::IntegerMatrix& lib,
                                const Rcpp::IntegerMatrix& pred,
                                const Rcpp::IntegerVector& E,
                                const Rcpp::IntegerVector& b,
                                const Rcpp::IntegerVector& tau,
                                int style = 1,
                                int zero_tolerance = 0,
                                int dist_metric = 2,
                                bool relative = true,
                                bool weighted = true,
                                int threads = 8,
                                int parallel_level = 0,
                                const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // Convert Rcpp::NumericMatrix to std::vector<std::vector<double>>
  int numRows = target.nrow();
  int numCols = target.ncol();
  std::vector<std::vector<double>> sourceMat(numRows, std::vector<double>(numCols));
  std::vector<std::vector<double>> targetMat(numRows, std::vector<double>(numCols));

  for (int r = 0; r < numRows; ++r) {
    for (int c = 0; c < numCols; ++c) {
      sourceMat[r][c] = source(r, c);
      targetMat[r][c] = target(r, c);
    }
  }

  // Initialize lib_indices and pred_indices
  std::vector<size_t> lib_indices;
  lib_indices.reserve(lib.nrow());
  std::vector<size_t> pred_indices;
  pred_indices.reserve(pred.nrow());

  // Convert lib and pred (1-based in R) to 0-based indices and set corresponding positions to true
  int currow;
  int curcol;
  int lib_col = lib.ncol();
  int pred_col = pred.ncol();

  if (lib_col == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
          lib_indices.push_back(static_cast<size_t>(lib(i,0)-1));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      // Convert to 0-based index
      currow = lib(i,0);
      curcol = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
          lib_indices.push_back(static_cast<size_t>(LocateGridIndices(currow, curcol, numRows, numCols)));
      }
    }
  }

  if (pred_col == 1) {
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(targetMat[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
          pred_indices.push_back(static_cast<size_t>(pred(i,0)-1));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      // Convert to 0-based index
      currow = pred(i,0);
      curcol = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(sourceMat[currow-1][curcol-1]) &&
          !std::isnan(targetMat[currow-1][curcol-1])){
          pred_indices.push_back(static_cast<size_t>(LocateGridIndices(currow, curcol, numRows, numCols)));
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  // Check the validity of the neighbor numbers
  std::vector<int> b_std;
  b_std.reserve(b.size());
  for (int i = 0; i < b.size(); ++i){
    if (b[i] > static_cast<int>(lib_indices.size())) {
      Rcpp::stop("Neighbor numbers count out of acceptable range at position %d (value: %d)", i + 1, b[i]);
    }
    b_std.push_back(b[i]);
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  std::vector<std::vector<double>> res_std = PC4Grid(
    sourceMat,
    targetMat,
    lib_indices,
    pred_indices,
    E_std,
    b_std,
    tau_std,
    style,
    zero_tolerance,
    dist_metric,
    relative,
    weighted,
    threads,
    parallel_level,
    dir_std);

  size_t n_rows = res_std.size();
  size_t n_cols = res_std[0].size();

  // Create an Rcpp::NumericMatrix with the same dimensions
  Rcpp::NumericMatrix result(n_rows, n_cols);

  // Fill the Rcpp::NumericMatrix with data from res_std
  for (size_t i = 0; i < n_rows; ++i) {
    for (size_t j = 0; j < n_cols; ++j) {
      result(i, j) = res_std[i][j];
    }
  }

  // Set column names for the result matrix
  Rcpp::colnames(result) = Rcpp::CharacterVector::create("E", "k", "tau", "positive", "negative", "dark");
  return result;
}

// Wrapper function to perform GCCM for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppGCCM4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    int E = 3,
    int tau = 1,
    int b = 5,
    bool simplex = true,
    double theta = 0,
    int threads = 8,
    int parallel_level = 0,
    int style = 1,
    int stack = 0,
    int dist_metric = 2,
    bool dist_average = true,
    bool single_sig = true,
    bool progressbar = false,
    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0),
    const Rcpp::NumericVector& win_ratio = Rcpp::NumericVector::create(0,0)) {
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  for (int i = 0; i < xMatrix.nrow(); ++i) {
    for (int j = 0; j < xMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert libsizes to a fundamental C++ data type
  int libsizes_dim = libsizes.ncol();
  std::vector<int> libsizes_cpp1;
  libsizes_cpp1.reserve(libsizes.nrow());
  std::vector<std::vector<int>> libsizes_cpp2(2);
  libsizes_cpp2[0].reserve(libsizes.nrow());
  libsizes_cpp2[1].reserve(libsizes.nrow());
  if (libsizes_dim == 1){
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_cpp1.push_back(libsizes(i, 0));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_cpp2[0].push_back(libsizes(i, 0));
      libsizes_cpp2[1].push_back(libsizes(i, 1));
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<int> lib_cpp1;
  lib_cpp1.reserve(lib.nrow());
  std::vector<std::pair<int, int>> lib_cpp2;
  lib_cpp2.reserve(lib.nrow());
  if (libsizes_dim == 1){
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
          lib_cpp1.push_back(lib(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // Convert to 0-based index
        int currow = lib(i,0);
        int curcol = lib(i,1);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          lib_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          lib_cpp2.emplace_back(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1])){
          lib_cpp2.emplace_back(lib(i, 0) - 1, lib(i, 1) - 1);
        }
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<int> pred_cpp1;
  pred_cpp1.reserve(pred.nrow());
  std::vector<std::pair<int, int>> pred_cpp2;
  pred_cpp2.reserve(pred.nrow());
  if (libsizes_dim == 1){
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
          pred_cpp1.push_back(pred(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // Convert to 0-based index
        int currow = pred(i,0);
        int curcol = pred(i,1);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          pred_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          pred_cpp2.emplace_back(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1])){
          pred_cpp2.emplace_back(pred(i, 0) - 1, pred(i, 1) - 1);
        }
      }
    }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0–8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_cpp = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_cpp.begin(), dir_cpp.end());
  dir_cpp.erase(std::unique(dir_cpp.begin(), dir_cpp.end()), dir_cpp.end());

  // if 0 is present, override all others
  if (std::find(dir_cpp.begin(), dir_cpp.end(), 0) != dir_cpp.end()) {
    dir_cpp = {0};
  }

  std::vector<double> win_ratio_cpp = Rcpp::as<std::vector<double>>(win_ratio);

  // Call the C++ function GCCM4Grid or GCCM4GridOneDim
  std::vector<std::vector<double>> result;
  if (libsizes_dim == 1){
    result = GCCM4GridOneDim(
      xMatrix_cpp,
      yMatrix_cpp,
      libsizes_cpp1,
      lib_cpp1,
      pred_cpp1,
      E,
      tau,
      b,
      simplex,
      theta,
      threads,
      parallel_level,
      style,
      stack,
      dist_metric,
      dist_average,
      single_sig,
      progressbar,
      dir_cpp
    );
  } else{
    result = GCCM4Grid(
      xMatrix_cpp,
      yMatrix_cpp,
      libsizes_cpp2,
      lib_cpp2,
      pred_cpp2,
      E,
      tau,
      b,
      simplex,
      theta,
      threads,
      parallel_level,
      style,
      stack,
      dist_metric,
      dist_average,
      single_sig,
      progressbar,
      dir_cpp,
      win_ratio_cpp
    );
  }

  Rcpp::NumericMatrix resultMatrix(result.size(), 5);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create("libsizes",
                 "x_xmap_y_mean","x_xmap_y_sig",
                 "x_xmap_y_lower","x_xmap_y_upper");
  return resultMatrix;
}

// Wrapper function to perform SCPCM for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericMatrix RcppSCPCM4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::NumericMatrix& zMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    const Rcpp::IntegerVector& b,
    bool simplex = true,
    double theta = 0,
    int threads = 8,
    int parallel_level = 0,
    bool cumulate = false,
    int style = 1,
    int stack = 0,
    int dist_metric = 2,
    bool dist_average = true,
    bool single_sig = true,
    bool progressbar = false,
    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0),
    const Rcpp::NumericVector& win_ratio = Rcpp::NumericVector::create(0,0)) {
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  for (int i = 0; i < xMatrix.nrow(); ++i) {
    for (int j = 0; j < xMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector of std::vectors
  std::vector<std::vector<double>> zMatrix_cpp(zMatrix.ncol());
  for (int i = 0; i < zMatrix.ncol(); ++i) {
    Rcpp::NumericVector covvar = zMatrix.column(i);
    zMatrix_cpp[i] = Rcpp::as<std::vector<double>>(covvar);
  }

  // Convert libsizes to a fundamental C++ data type
  int libsizes_dim = libsizes.ncol();
  std::vector<int> libsizes_cpp1;
  libsizes_cpp1.reserve(libsizes.nrow());
  std::vector<std::vector<int>> libsizes_cpp2(2);
  libsizes_cpp2[0].reserve(libsizes.nrow());
  libsizes_cpp2[1].reserve(libsizes.nrow());
  if (libsizes_dim == 1){
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_cpp1.push_back(libsizes(i, 0));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_cpp2[0].push_back(libsizes(i, 0));
      libsizes_cpp2[1].push_back(libsizes(i, 1));
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<int> lib_cpp1;
  lib_cpp1.reserve(lib.nrow());
  std::vector<std::pair<int, int>> lib_cpp2;
  lib_cpp2.reserve(lib.nrow());
  if (libsizes_dim == 1){
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(lib(i,0)-1) / numCols][(lib(i,0)-1) % numCols])){
          lib_cpp1.push_back(lib(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // Convert to 0-based index
        int currow = lib(i,0);
        int curcol = lib(i,1);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          lib_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (lib_dim == 1){
      for (int i = 0; i < lib.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          lib_cpp2.emplace_back(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < lib.nrow(); ++i) {
        // disallow lib indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[lib(i, 0) - 1][lib(i, 1) - 1])){
          lib_cpp2.emplace_back(lib(i, 0) - 1, lib(i, 1) - 1);
        }
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<int> pred_cpp1;
  pred_cpp1.reserve(pred.nrow());
  std::vector<std::pair<int, int>> pred_cpp2;
  pred_cpp2.reserve(pred.nrow());
  if (libsizes_dim == 1){
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
            !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
          pred_cpp1.push_back(pred(i, 0) - 1);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // Convert to 0-based index
        int currow = pred(i,0);
        int curcol = pred(i,1);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[currow-1][curcol-1]) &&
            !std::isnan(yMatrix_cpp[currow-1][curcol-1])){
          pred_cpp1.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
        }
      }
    }
  } else {
    if (pred_dim == 1){
      for (int i = 0; i < pred.nrow(); ++i) {
        std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[rowcolnum[0]][rowcolnum[1]]) &&
            !std::isnan(yMatrix_cpp[rowcolnum[0]][rowcolnum[1]])){
          pred_cpp2.emplace_back(rowcolnum[0], rowcolnum[1]);
        }
      }
    } else {
      for (int i = 0; i < pred.nrow(); ++i) {
        // disallow pred indices to point to vectors with NaN
        if (!std::isnan(xMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1]) &&
            !std::isnan(yMatrix_cpp[pred(i, 0) - 1][pred(i, 1) - 1])){
          pred_cpp2.emplace_back(pred(i, 0) - 1, pred(i, 1) - 1);
        }
      }
    }
  }

  // Convert Rcpp::IntegerVector to std::vector<int>
  std::vector<int> E_cpp = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_cpp = Rcpp::as<std::vector<int>>(tau);
  std::vector<int> b_cpp = Rcpp::as<std::vector<int>>(b);

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_cpp = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_cpp.begin(), dir_cpp.end());
  dir_cpp.erase(std::unique(dir_cpp.begin(), dir_cpp.end()), dir_cpp.end());

  // if 0 is present, override all others
  if (std::find(dir_cpp.begin(), dir_cpp.end(), 0) != dir_cpp.end()) {
    dir_cpp = {0};
  }

  std::vector<double> win_ratio_cpp = Rcpp::as<std::vector<double>>(win_ratio);

  // Call the C++ function SCPCM4Grid or SCPCM4GridOneDim
  std::vector<std::vector<double>> result;
  if (libsizes_dim == 1){
    result = SCPCM4GridOneDim(
      xMatrix_cpp,
      yMatrix_cpp,
      zMatrix_cpp,
      libsizes_cpp1,
      lib_cpp1,
      pred_cpp1,
      E_cpp,
      tau_cpp,
      b_cpp,
      simplex,
      theta,
      threads,
      parallel_level,
      cumulate,
      style,
      stack,
      dist_metric,
      dist_average,
      single_sig,
      progressbar,
      dir_cpp
    );
  } else{
    result = SCPCM4Grid(
      xMatrix_cpp,
      yMatrix_cpp,
      zMatrix_cpp,
      libsizes_cpp2,
      lib_cpp2,
      pred_cpp2,
      E_cpp,
      tau_cpp,
      b_cpp,
      simplex,
      theta,
      threads,
      parallel_level,
      cumulate,
      style,
      stack,
      dist_metric,
      dist_average,
      single_sig,
      progressbar,
      dir_cpp,
      win_ratio_cpp
    );
  }

  // Convert std::vector<std::vector<double>> to Rcpp::NumericMatrix
  Rcpp::NumericMatrix resultMatrix(result.size(), 9);
  for (size_t i = 0; i < result.size(); ++i) {
    resultMatrix(i, 0) = result[i][0];
    resultMatrix(i, 1) = result[i][1];
    resultMatrix(i, 2) = result[i][2];
    resultMatrix(i, 3) = result[i][3];
    resultMatrix(i, 4) = result[i][4];
    resultMatrix(i, 5) = result[i][5];
    resultMatrix(i, 6) = result[i][6];
    resultMatrix(i, 7) = result[i][7];
    resultMatrix(i, 8) = result[i][8];
  }

  // Set column names for the result matrix
  Rcpp::colnames(resultMatrix) = Rcpp::CharacterVector::create(
    "libsizes","T_mean","D_mean",
    "T_sig","T_lower","T_upper",
    "D_sig","D_lower","D_upper");
  return resultMatrix;
}

// Wrapper function to perform GCMC for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGCMC4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int b = 4,
    int r = 0,
    int style = 1,
    int dist_metric = 2,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = false,
    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)){
  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  for (int i = 0; i < xMatrix.nrow(); ++i) {
    for (int j = 0; j < xMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
    }
  }

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  size_t validCellNum = 0;
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      if (!std::isnan(yMatrix(i, j))){
        validCellNum += 1;
      }
      yMatrix_cpp[i][j] = yMatrix(i, j);
    }
  }

  // Convert Rcpp IntegerVector to std::vector<int>
  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  // Convert Rcpp IntegerMatrix to std::vector<int>
  int n_libcol = lib.ncol();
  int n_predcol = pred.ncol();
  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();

  std::vector<size_t> lib_std;
  lib_std.reserve(lib.nrow());
  if (n_libcol == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        lib_std.push_back(static_cast<size_t>(lib(i,0) - 1));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int rowLibIndice = lib(i,0);
      int colLibIndice = lib(i,1);
      // disallow lib indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[rowLibIndice-1][colLibIndice-1]) &&
          !std::isnan(yMatrix_cpp[rowLibIndice-1][colLibIndice-1])){
        lib_std.push_back(static_cast<size_t>(LocateGridIndices(rowLibIndice, colLibIndice, numRows, numCols)));
      }
    }
  }

  std::vector<size_t> pred_std;
  pred_std.reserve(pred.nrow());
  if (n_predcol == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols]) &&
          !std::isnan(yMatrix_cpp[(pred(i,0)-1) / numCols][(pred(i,0)-1) % numCols])){
        pred_std.push_back(static_cast<size_t>(pred(i, 0) - 1));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int rowPredIndice = pred(i,0);
      int colPredIndice = pred(i,1);
      // disallow pred indices to point to vectors with NaN
      if (!std::isnan(xMatrix_cpp[rowPredIndice-1][colPredIndice-1]) &&
          !std::isnan(yMatrix_cpp[rowPredIndice-1][colPredIndice-1])){
        pred_std.push_back(static_cast<size_t>(LocateGridIndices(rowPredIndice, colPredIndice, numRows, numCols)));
      }
    }
  }

  if (b <= 3 || static_cast<size_t>(b) > validCellNum) {
    Rcpp::stop("k must be greater than 3 and no larger than the number of non-NA values.\n"
               "An empirical rule of thumb is to set k = sqrt(E * N),\n"
               "where E is the embedding dimension and N is the number of valid observations in the prediction set.");
  }

  // -- Convert libsizes --------------------------------------------------

  int libsizes_dim = libsizes.ncol();
  std::vector<size_t> libsizes_std;
  libsizes_std.reserve(libsizes.nrow());
  if (libsizes_dim == 1){
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_std.push_back(static_cast<size_t>(libsizes(i, 0)));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_std.push_back(static_cast<size_t>(libsizes(i, 0) * libsizes(i, 1)));
    }
  }

  for (size_t& libsize : libsizes_std){
    if (libsize < static_cast<size_t>(b+r) || libsize > lib_std.size()) libsize = lib_std.size();
  }
  std::sort(libsizes_std.begin(), libsizes_std.end());
  libsizes_std.erase(std::unique(libsizes_std.begin(), libsizes_std.end()), libsizes_std.end());

  if (libsizes_std.empty()) {
    Rcpp::warning("[Warning] No valid libsizes after filtering. Using full library size as fallback.");
    libsizes_std.push_back(lib_std.size());
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_cpp = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_cpp.begin(), dir_cpp.end());
  dir_cpp.erase(std::unique(dir_cpp.begin(), dir_cpp.end()), dir_cpp.end());

  // if 0 is present, override all others
  if (std::find(dir_cpp.begin(), dir_cpp.end(), 0) != dir_cpp.end()) {
    dir_cpp = {0};
  }

  // Generate embeddings
  std::vector<std::vector<double>> e1 = GenGridEmbeddings(xMatrix_cpp, E[0], tau_std[0], style, dir_cpp);
  std::vector<std::vector<double>> e2 = GenGridEmbeddings(yMatrix_cpp, E[1], tau_std[1], style, dir_cpp);

  // Perform GCMC for spatial grid data
  CMCRes res = CrossMappingCardinality(e1,e2,libsizes_std,lib_std,pred_std,
                                       static_cast<size_t>(b),static_cast<size_t>(r),
                                       dist_metric,threads,parallel_level,progressbar);

  // Convert mean_aucs to Rcpp::DataFrame
  std::vector<double> libs, aucs;
  for (const auto& cm : res.causal_strength) {
    libs.push_back(cm[0]);
    aucs.push_back(cm[1]);
  }

  Rcpp::DataFrame cs_df = Rcpp::DataFrame::create(
    Rcpp::Named("libsizes") = libs,
    Rcpp::Named("x_xmap_y_mean") = aucs
  );

  // Wrap causal_strength with names
  Rcpp::DataFrame xmap_df = Rcpp::DataFrame::create(
    Rcpp::Named("neighbors") = res.cross_mapping[0],
    Rcpp::Named("x_xmap_y_mean") = res.cross_mapping[1],
    Rcpp::Named("x_xmap_y_sig") = res.cross_mapping[2],
    Rcpp::Named("x_xmap_y_lower") = res.cross_mapping[3],
    Rcpp::Named("x_xmap_y_upper")  = res.cross_mapping[4]
  );

  return Rcpp::List::create(
    Rcpp::Named("xmap") = xmap_df,
    Rcpp::Named("cs") = cs_df
  );
}

// Wrapper function to perform Geographical Pattern Causality (GPC) for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::List RcppGPC4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int style = 1,
    int b = 4,
    int zero_tolerance = 0,
    int dist_metric = 2,
    bool relative = true,
    bool weighted = true,
    int threads = 8,
    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // --- Convert inputs to C++ types ---
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));

  size_t validCellNum = 0;
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
      yMatrix_cpp[i][j] = yMatrix(i, j);
      if (!std::isnan(yMatrix(i, j))) validCellNum += 1;
    }
  }

  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();
  int n_libcol = lib.ncol();
  int n_predcol = pred.ncol();

  // --- Convert library and prediction indices ---
  std::vector<size_t> lib_std;
  lib_std.reserve(lib.nrow());

  if (n_libcol == 1) {
    for (int i = 0; i < lib.nrow(); ++i) {
      int idx = lib(i, 0) - 1;
      int r = idx / numCols, c = idx % numCols;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        lib_std.push_back(static_cast<size_t>(idx));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int r = lib(i, 0) - 1, c = lib(i, 1) - 1;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        lib_std.push_back(static_cast<size_t>(LocateGridIndices(r + 1, c + 1, numRows, numCols)));
      }
    }
  }

  std::vector<size_t> pred_std;
  pred_std.reserve(pred.nrow());

  if (n_predcol == 1) {
    for (int i = 0; i < pred.nrow(); ++i) {
      int idx = pred(i, 0) - 1;
      int r = idx / numCols, c = idx % numCols;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        pred_std.push_back(static_cast<size_t>(idx));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int r = pred(i, 0) - 1, c = pred(i, 1) - 1;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        pred_std.push_back(static_cast<size_t>(LocateGridIndices(r + 1, c + 1, numRows, numCols)));
      }
    }
  }

  // Sort + unique lib/pred
  std::sort(lib_std.begin(), lib_std.end());
  lib_std.erase(
    std::unique(lib_std.begin(), lib_std.end()),
    lib_std.end()
  );

  std::sort(pred_std.begin(), pred_std.end());
    pred_std.erase(
      std::unique(pred_std.begin(), pred_std.end()),
      pred_std.end()
  );

  // Copy prediction indices for mapping back to original dataset
  std::vector<size_t> pred_indices = pred_std;

  // --- Validate parameters ---
  if (b < 2 || static_cast<size_t>(b) > validCellNum)
    Rcpp::stop("k cannot be less than or equal to 2 or greater than the number of non-NA values.");

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // --- Generate embeddings ---
  std::vector<std::vector<double>> Mx = GenGridEmbeddings(xMatrix_cpp, E_std[0], tau_std[0], style, dir_std);
  std::vector<std::vector<double>> My = GenGridEmbeddings(yMatrix_cpp, E_std[1], tau_std[1], style, dir_std);

  // --- Prepare for data slicing ---
  std::vector<size_t> selected_indices;
  selected_indices.reserve(lib_std.size() + pred_std.size());
  for (size_t i = 0; i < lib_std.size(); ++i)
    selected_indices.push_back(lib_std[i]);
  for (size_t i = 0; i < pred_std.size(); ++i)
    selected_indices.push_back(pred_std[i]);
  std::sort(selected_indices.begin(), selected_indices.end());
  selected_indices.erase(
    std::unique(selected_indices.begin(), selected_indices.end()),
    selected_indices.end()
  );

  // --- Check if full set is used ---
  bool use_subset = (selected_indices.size() < Mx.size());

  // --- Perform GPC analysis ---
  PatternCausalityRes res;

  if (!use_subset) {
    // --- Full data: no slicing needed ---
    res = PatternCausality(
      Mx, My, lib_std, pred_std, b, zero_tolerance,
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

    // --- Subset mode: build index map ---
    std::unordered_map<size_t, size_t> index_map;
    index_map.reserve(selected_indices.size());

    for (size_t i = 0; i < selected_indices.size(); ++i) {
      index_map[selected_indices[i]] = i;
    }

    // --- Remap lib indices ---
    for (size_t i = 0; i < lib_std.size(); ++i) {
      lib_std[i] = index_map[lib_std[i]];
    }

    // --- Remap pred indices ---
    for (size_t i = 0; i < pred_std.size(); ++i) {
      pred_std[i] = index_map[pred_std[i]];
    }

    // --- Run pattern causality on subset ---
    res = PatternCausality(
      Mx_sub, My_sub, lib_std, pred_std, b, zero_tolerance,
      dist_metric, relative, weighted, threads);
  }

  // --- Convert pattern matrix to Rcpp::NumericMatrix ---
  size_t nrow = res.matrice.size();
  size_t ncol = nrow > 0 ? res.matrice[0].size() : 0;
  Rcpp::NumericMatrix pattern_mat(nrow, ncol);

  for (size_t i = 0; i < nrow; ++i) {
    for (size_t j = 0; j < ncol; ++j) {
      pattern_mat(i, j) = res.matrice[i][j];
    }
  }
  
  // Assign row and column names if available
  if (res.PatternStrings.size() == ncol) {
    pattern_mat.attr("dimnames") = Rcpp::List::create(
      Rcpp::CharacterVector(res.PatternStrings.begin(), res.PatternStrings.end()),
      Rcpp::CharacterVector(res.PatternStrings.begin(), res.PatternStrings.end())
    );
  }

  // --- Create DataFrame for per-sample causality ---
  size_t n_samples = res.NoCausality.size();
  Rcpp::IntegerVector real_index(n_samples);  // original indices (+1 for R)
  Rcpp::CharacterVector pattern_labels(n_samples, "no");

  for (size_t rl = 0; rl < n_samples; ++rl) {
    // Restore original index (+1 for R)
    real_index[rl] = static_cast<int>(pred_indices[rl] + 1);

    // Map pattern_types (0–3) → descriptive string labels
    switch (res.PatternTypes[rl]) {
      case 0: pattern_labels[rl]  = "no"; break;
      case 1: pattern_labels[rl]  = "positive"; break;
      case 2: pattern_labels[rl]  = "negative"; break;
      case 3: pattern_labels[rl]  = "dark"; break;
      default: pattern_labels[rl] = "unknown"; break;
    }
  }

  // --- Build DataFrame outputs ---
  Rcpp::DataFrame causality_df = Rcpp::DataFrame::create(
    Rcpp::Named("index") = real_index,
    Rcpp::Named("no") = Rcpp::NumericVector(res.NoCausality.begin(), res.NoCausality.end()),
    Rcpp::Named("positive") = Rcpp::NumericVector(res.PositiveCausality.begin(), res.PositiveCausality.end()),
    Rcpp::Named("negative") = Rcpp::NumericVector(res.NegativeCausality.begin(), res.NegativeCausality.end()),
    Rcpp::Named("dark") = Rcpp::NumericVector(res.DarkCausality.begin(), res.DarkCausality.end()),
    Rcpp::Named("type") = pattern_labels
  );

  Rcpp::DataFrame summary_df = Rcpp::DataFrame::create(
    Rcpp::Named("type") = Rcpp::CharacterVector::create("positive", "negative", "dark"),
    Rcpp::Named("strength") = Rcpp::NumericVector::create(res.TotalPos, res.TotalNeg, res.TotalDark)
  );

  // --- Return structured result ---
  return Rcpp::List::create(
    Rcpp::Named("causality") = causality_df,
    Rcpp::Named("summary") = summary_df,
    Rcpp::Named("pattern") = pattern_mat
  );
}

// Wrapper function to perform Robust Geographical Pattern Causality for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::DataFrame RcppGPCRobust4Grid(
    const Rcpp::NumericMatrix& xMatrix,
    const Rcpp::NumericMatrix& yMatrix,
    const Rcpp::IntegerMatrix& libsizes,
    const Rcpp::IntegerMatrix& lib,
    const Rcpp::IntegerMatrix& pred,
    const Rcpp::IntegerVector& E,
    const Rcpp::IntegerVector& tau,
    int style = 1,
    int b = 4,
    int boot = 99,
    bool replace_sampling = false,
    unsigned long long seed = 42,
    int zero_tolerance = 0,
    int dist_metric = 2,
    bool relative = true,
    bool weighted = true,
    int threads = 8,
    int parallel_level = 0,
    bool progressbar = false,
    const Rcpp::IntegerVector& dir = Rcpp::IntegerVector::create(0)) {
  // --- Convert inputs to C++ types ---
  std::vector<std::vector<double>> xMatrix_cpp(xMatrix.nrow(), std::vector<double>(xMatrix.ncol()));
  std::vector<std::vector<double>> yMatrix_cpp(yMatrix.nrow(), std::vector<double>(yMatrix.ncol()));

  size_t validCellNum = 0;
  for (int i = 0; i < yMatrix.nrow(); ++i) {
    for (int j = 0; j < yMatrix.ncol(); ++j) {
      xMatrix_cpp[i][j] = xMatrix(i, j);
      yMatrix_cpp[i][j] = yMatrix(i, j);
      if (!std::isnan(yMatrix(i, j))) validCellNum += 1;
    }
  }

  std::vector<int> E_std = Rcpp::as<std::vector<int>>(E);
  std::vector<int> tau_std = Rcpp::as<std::vector<int>>(tau);

  int numRows = yMatrix.nrow();
  int numCols = yMatrix.ncol();
  int n_libcol = lib.ncol();
  int n_predcol = pred.ncol();

  // --- Convert library and prediction indices ---
  std::vector<size_t> lib_std;
  lib_std.reserve(lib.nrow());

  if (n_libcol == 1) {
    for (int i = 0; i < lib.nrow(); ++i) {
      int idx = lib(i, 0) - 1;
      int r = idx / numCols, c = idx % numCols;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        lib_std.push_back(static_cast<size_t>(idx));
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      int r = lib(i, 0) - 1, c = lib(i, 1) - 1;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        lib_std.push_back(static_cast<size_t>(LocateGridIndices(r + 1, c + 1, numRows, numCols)));
      }
    }
  }

  std::vector<size_t> pred_std;
  pred_std.reserve(pred.nrow());

  if (n_predcol == 1) {
    for (int i = 0; i < pred.nrow(); ++i) {
      int idx = pred(i, 0) - 1;
      int r = idx / numCols, c = idx % numCols;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        pred_std.push_back(static_cast<size_t>(idx));
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      int r = pred(i, 0) - 1, c = pred(i, 1) - 1;
      if (!std::isnan(xMatrix_cpp[r][c]) && !std::isnan(yMatrix_cpp[r][c])) {
        pred_std.push_back(static_cast<size_t>(LocateGridIndices(r + 1, c + 1, numRows, numCols)));
      }
    }
  }

  // Sort + unique lib/pred
  std::sort(lib_std.begin(), lib_std.end());
  lib_std.erase(
    std::unique(lib_std.begin(), lib_std.end()),
    lib_std.end()
  );

  std::sort(pred_std.begin(), pred_std.end());
    pred_std.erase(
      std::unique(pred_std.begin(), pred_std.end()),
      pred_std.end()
  );

  // -- Convert libsizes ---
  int libsizes_dim = libsizes.ncol();
  std::vector<size_t> libsizes_std;
  libsizes_std.reserve(libsizes.nrow());
  if (libsizes_dim == 1) {
    for (int i = 0; i < libsizes.nrow(); ++i) {
      libsizes_std.push_back(static_cast<size_t>(libsizes(i, 0)));
    }
  } else {
    for (int i = 0; i < libsizes.nrow(); ++i) { // Fill all the sub-vector
      libsizes_std.push_back(static_cast<size_t>(libsizes(i, 0) * libsizes(i, 1)));
    }
  }

  std::vector<size_t> valid_libsizes;
  valid_libsizes.reserve(libsizes_std.size());
  for (size_t s : libsizes_std) {
    if (s > static_cast<size_t>(b) && s <= lib_std.size())
      valid_libsizes.push_back(s);
  }

  std::sort(valid_libsizes.begin(), valid_libsizes.end());
  valid_libsizes.erase(std::unique(valid_libsizes.begin(), valid_libsizes.end()), valid_libsizes.end());

  if (valid_libsizes.empty()) {
    Rcpp::warning("[Warning] No valid libsizes after filtering. Using full library size as fallback.");
    valid_libsizes.push_back(lib_std.size());
  }

  // --- Validate parameters ---
  if (b < 2 || static_cast<size_t>(b) > validCellNum)
    Rcpp::stop("k cannot be less than or equal to 2 or greater than the number of non-NA values.");

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_cpp = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_cpp.begin(), dir_cpp.end());
  dir_cpp.erase(std::unique(dir_cpp.begin(), dir_cpp.end()), dir_cpp.end());

  // if 0 is present, override all others
  if (std::find(dir_cpp.begin(), dir_cpp.end(), 0) != dir_cpp.end()) {
    dir_cpp = {0};
  }

  // --- Generate embeddings ---
  std::vector<std::vector<double>> Mx = GenGridEmbeddings(xMatrix_cpp, E_std[0], tau_std[0], style, dir_cpp);
  std::vector<std::vector<double>> My = GenGridEmbeddings(yMatrix_cpp, E_std[1], tau_std[1], style, dir_cpp);

  // --- Prepare for data slicing ---
  std::vector<size_t> selected_indices;
  selected_indices.reserve(lib_std.size() + pred_std.size());
  for (size_t i = 0; i < lib_std.size(); ++i)
    selected_indices.push_back(lib_std[i]);
  for (size_t i = 0; i < pred_std.size(); ++i)
    selected_indices.push_back(pred_std[i]);
  std::sort(selected_indices.begin(), selected_indices.end());
  selected_indices.erase(
    std::unique(selected_indices.begin(), selected_indices.end()),
    selected_indices.end()
  );

  // --- Check if full set is used ---
  bool use_subset = (selected_indices.size() < Mx.size());

  // --- Perform Robust GPC analysis ---
  std::vector<std::vector<std::vector<double>>> res;

  if (!use_subset) {
    // --- Full data: no slicing needed ---
    res = RobustPatternCausality(
      Mx, My, valid_libsizes, lib_std, pred_std, b, boot, replace_sampling, seed, zero_tolerance,
      dist_metric, relative, weighted, threads, parallel_level, progressbar);
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

    // --- Subset mode: build index map ---
    std::unordered_map<size_t, size_t> index_map;
    index_map.reserve(selected_indices.size());

    for (size_t i = 0; i < selected_indices.size(); ++i) {
      index_map[selected_indices[i]] = i;
    }

    // --- Remap lib indices ---
    for (size_t i = 0; i < lib_std.size(); ++i) {
      lib_std[i] = index_map[lib_std[i]];
    }

    // --- Remap pred indices ---
    for (size_t i = 0; i < pred_std.size(); ++i) {
      pred_std[i] = index_map[pred_std[i]];
    }

    // --- Run pattern causality on subset ---
    res = RobustPatternCausality(
      Mx_sub, My_sub, valid_libsizes, lib_std, pred_std, b, boot, replace_sampling, seed, zero_tolerance,
      dist_metric, relative, weighted, threads, parallel_level, progressbar);
  }

  // --- Result Processing ---

  // res structure: [3][libsizes][boot]
  // dimension 0: metric type (0=Positive,1=Negative,2=Dark)
  // dimension 1: libsizes index
  // dimension 2: bootstrap replicates

  int n_types = 3;
  int n_libsizes = static_cast<int>(valid_libsizes.size());
  int n_boot = static_cast<int>(res[0][0].size());

  // Prepare vectors to hold dataframe columns
  std::vector<size_t> df_libsizes;
  std::vector<std::string> df_type;
  std::vector<double> df_causality;
  std::vector<double> df_q05, df_q50, df_q95;  // For quantiles if boot > 1

  bool has_bootstrap = (n_boot > 1);
  const std::string types[3] = {"positive", "negative", "dark"};

  if (!has_bootstrap) {
    // boot == 1, simple long format: columns = libsizes, type, causality
    df_libsizes.reserve(n_types * n_libsizes);
    df_type.reserve(n_types * n_libsizes);
    df_causality.reserve(n_types * n_libsizes);

    for (int t = 0; t < n_types; ++t) {
      for (int l = 0; l < n_libsizes; ++l) {
        df_libsizes.push_back(valid_libsizes[l]);
        df_type.push_back(types[t]);
        if(std::isnan(res[t][l][0])) res[t][l][0] = 0; // replace nan causal strength with 0
        df_causality.push_back(res[t][l][0]);
      }
    }

    return Rcpp::DataFrame::create(
      Rcpp::Named("libsizes") = df_libsizes,
      Rcpp::Named("type") = df_type,
      Rcpp::Named("causality") = df_causality
    );
  } else {
    // boot > 1, summary with mean and quantiles
    df_libsizes.reserve(n_types * n_libsizes);
    df_type.reserve(n_types * n_libsizes);
    df_causality.reserve(n_types * n_libsizes);
    df_q05.reserve(n_types * n_libsizes);
    df_q50.reserve(n_types * n_libsizes);
    df_q95.reserve(n_types * n_libsizes);

    for (int t = 0; t < n_types; ++t) {
      for (int l = 0; l < n_libsizes; ++l) {
        const std::vector<double>& boot_vals = res[t][l];
        double mean_val = CppMean(boot_vals, true);
        std::vector<double> qs = CppQuantile(boot_vals, {0.05, 0.5, 0.95}, true);

        // replace nan causal strength with 0
        if(std::isnan(mean_val)) mean_val = 0;
        for(double& q : qs){
          if(std::isnan(q)) q = 0;
        }

        df_libsizes.push_back(valid_libsizes[l]);
        df_type.push_back(types[t]);
        df_causality.push_back(mean_val);
        df_q05.push_back(qs[0]);
        df_q50.push_back(qs[1]);
        df_q95.push_back(qs[2]);
      }
    }

    return Rcpp::DataFrame::create(
      Rcpp::Named("libsizes") = df_libsizes,
      Rcpp::Named("type") = df_type,
      Rcpp::Named("mean") = df_causality,
      Rcpp::Named("q05") = df_q05,
      Rcpp::Named("q50") = df_q50,
      Rcpp::Named("q95") = df_q95
    );
  }
}

// Wrapper function to perform SGC for spatial grid data without bootstrapped significance
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSGCSingle4Grid(const Rcpp::NumericMatrix& x,
                                       const Rcpp::NumericMatrix& y,
                                       const Rcpp::IntegerVector& dir,
                                       const Rcpp::IntegerMatrix& lib,
                                       const Rcpp::IntegerMatrix& pred,
                                       int k,
                                       double base = 2,
                                       bool symbolize = true,
                                       bool normalize = false){
  int numRows = y.nrow();
  int numCols = y.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xmat(x.nrow(), std::vector<double>(x.ncol()));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      xmat[i][j] = x(i, j);
    }
  }
  std::vector<std::vector<double>> ymat(y.nrow(), std::vector<double>(y.ncol()));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      ymat[i][j] = y(i, j);
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<std::pair<int, int>> lib_std;
  lib_std.reserve(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        lib_std.emplace_back(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      if (!std::isnan(xmat[lib(i, 0) - 1][lib(i, 1) - 1]) &&
          !std::isnan(ymat[lib(i, 0) - 1][lib(i, 1) - 1])){
        lib_std.emplace_back(lib(i, 0) - 1, lib(i, 1) - 1);
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<std::pair<int, int>> pred_std;
  pred_std.reserve(pred.nrow());

  if (pred_dim == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        pred_std.emplace_back(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      if (!std::isnan(xmat[pred(i, 0) - 1][pred(i, 1) - 1]) &&
          !std::isnan(ymat[pred(i, 0) - 1][pred(i, 1) - 1])){
        pred_std.emplace_back(pred(i, 0) - 1, pred(i, 1) - 1);
      }
    }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // Perform SGC for spatial grid data
  std::vector<double> sc = SGCSingle4Grid(
    xmat,
    ymat,
    dir_std,
    lib_std,
    pred_std,
    k,
    base,
    symbolize,
    normalize
  );

  // Convert the result back to Rcpp::NumericVector
  Rcpp::NumericVector sc_res = Rcpp::wrap(sc);
  sc_res.names() = Rcpp::CharacterVector::create(
    "statistic for x → y causality",
    "statistic for y → x causality"
  );

  return sc_res;
}

// Wrapper function to perform SGC for spatial grid data
// [[Rcpp::export(rng = false)]]
Rcpp::NumericVector RcppSGC4Grid(const Rcpp::NumericMatrix& x,
                                 const Rcpp::NumericMatrix& y,
                                 const Rcpp::IntegerVector& dir,
                                 const Rcpp::IntegerMatrix& lib,
                                 const Rcpp::IntegerMatrix& pred,
                                 const Rcpp::IntegerMatrix& block,
                                 int k,
                                 int threads = 8,
                                 int boot = 399,
                                 double base = 2,
                                 unsigned long long seed = 42,
                                 bool symbolize = true,
                                 bool normalize = false,
                                 bool progressbar = true){
  int numRows = y.nrow();
  int numCols = y.ncol();

  // Convert Rcpp NumericMatrix to std::vector<std::vector<double>>
  std::vector<std::vector<double>> xmat(x.nrow(), std::vector<double>(x.ncol()));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      xmat[i][j] = x(i, j);
    }
  }
  std::vector<std::vector<double>> ymat(y.nrow(), std::vector<double>(y.ncol()));
  for (int i = 0; i < numRows; ++i) {
    for (int j = 0; j < numCols; ++j) {
      ymat[i][j] = y(i, j);
    }
  }

  // Convert lib to a fundamental C++ data type
  int lib_dim = lib.ncol();
  std::vector<std::pair<int, int>> lib_std;
  lib_std.reserve(lib.nrow());

  if (lib_dim == 1){
    for (int i = 0; i < lib.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(lib(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        lib_std.emplace_back(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < lib.nrow(); ++i) {
      if (!std::isnan(xmat[lib(i, 0) - 1][lib(i, 1) - 1]) &&
          !std::isnan(ymat[lib(i, 0) - 1][lib(i, 1) - 1])){
        lib_std.emplace_back(lib(i, 0) - 1, lib(i, 1) - 1);
      }
    }
  }

  // Convert pred to a fundamental C++ data type
  int pred_dim = pred.ncol();
  std::vector<std::pair<int, int>> pred_std;
  pred_std.reserve(pred.nrow());

  if (pred_dim == 1){
    for (int i = 0; i < pred.nrow(); ++i) {
      std::vector<int> rowcolnum = RowColFromGrid(pred(i, 0) - 1, numCols);
      if (!std::isnan(xmat[rowcolnum[0]][rowcolnum[1]]) &&
          !std::isnan(ymat[rowcolnum[0]][rowcolnum[1]])){
        pred_std.emplace_back(rowcolnum[0], rowcolnum[1]);
      }
    }
  } else {
    for (int i = 0; i < pred.nrow(); ++i) {
      if (!std::isnan(xmat[pred(i, 0) - 1][pred(i, 1) - 1]) &&
          !std::isnan(ymat[pred(i, 0) - 1][pred(i, 1) - 1])){
        pred_std.emplace_back(pred(i, 0) - 1, pred(i, 1) - 1);
      }
    }
  }

  // Convert block to a fundamental C++ data type
  int b_dim = block.ncol();
  std::vector<int> b_std;
  b_std.reserve(block.nrow());
  if (b_dim == 1){
    for (int i = 0; i < block.nrow(); ++i) {
      b_std.push_back(block(i, 0));
    }
  } else {
    for (int i = 0; i < block.nrow(); ++i) {
        // Convert to 0-based index
        int currow = block(i,0);
        int curcol = block(i,1);
        b_std.push_back(LocateGridIndices(currow, curcol, numRows, numCols));
      }
  }

  // check each element of dir before conversion
  for (int d : dir) {
    if (d < 0 || d > 8) {
      Rcpp::stop("direction vector elements must be in 0-8; 0 represents all directions, 1–8 correspond to NW,N,NE,W,E,SW,S,SE");
    }
  }

  // convert to std::vector<int>
  std::vector<int> dir_std = Rcpp::as<std::vector<int>>(dir);

  // remove duplicates
  std::sort(dir_std.begin(), dir_std.end());
  dir_std.erase(std::unique(dir_std.begin(), dir_std.end()), dir_std.end());

  // if 0 is present, override all others
  if (std::find(dir_std.begin(), dir_std.end(), 0) != dir_std.end()) {
    dir_std = {0};
  }

  // Perform SGC for spatial grid data
  std::vector<double> sc = SGC4Grid(
    xmat,
    ymat,
    dir_std,
    lib_std,
    pred_std,
    b_std,
    k,
    threads,
    boot,
    base,
    seed,
    symbolize,
    normalize,
    progressbar
  );

  // Convert the result back to Rcpp::NumericVector
  Rcpp::NumericVector sc_res = Rcpp::wrap(sc);
  sc_res.names() = Rcpp::CharacterVector::create(
    "statistic for x → y causality",
    "significance for x → y causality",
    "statistic for y → x causality",
    "significance for y → x causality"
  );

  return sc_res;
}

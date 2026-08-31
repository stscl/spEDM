#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <queue>
#include <unordered_set>
#include <set>
#include <limits>
#include <numeric>
#include <algorithm>
#include <utility>
#include <cstdint>
#include "CppStats.h"

/**
 * Converts a 2D grid position (row, column) to a 1D index in row-major order.
 * This function calculates the corresponding one-dimensional index based on the given
 * row and column position within a grid of specified dimensions.
 *
 * Parameters:
 *   curRow   - The current row number (1-based indexing).
 *   curCol   - The current column number (1-based indexing).
 *   totalRow - The total number of rows in the grid.
 *   totalCol - The total number of columns in the grid.
 *
 * Returns:
 *   The computed 1D index (0-based indexing) corresponding to the input row and column.
 */
int LocateGridIndices(int curRow, int curCol, int totalRow, int totalCol) {
  return (curRow - 1) * totalCol + curCol - 1;
}

/**
 * Converts a 1D grid index (cell number) to a 2D grid position (row, column) in row-major order.
 * This function determines the corresponding row and column indices based on a given
 * one-dimensional index within a grid of specified dimensions.
 *
 * Parameters:
 *   cellNum  - The 1D index of the cell (0-based indexing).
 *   totalCol - The total number of columns in the grid.
 *
 * Returns:
 *   A std::vector<int> containing the row index and column index (both 0-based indexing).
 */
std::vector<int> RowColFromGrid(int cellNum, int totalCol) {
  int row = cellNum / totalCol;
  int col = cellNum % totalCol;
  return {row, col};
}

/**
 * Converts a 2D grid data matrix (vector of vectors) to a 1D vector by concatenating the rows.
 * This function iterates over each row of the input matrix and appends the elements to the resulting vector.
 *
 * Parameters:
 *   Matrix - A 2D vector containing the grid data, where each element is a row of double values.
 *
 * Returns:
 *   A 1D vector containing all the elements of the input matrix, arranged row by row.
 */
std::vector<double> GridMat2Vec(const std::vector<std::vector<double>>& Matrix){
  std::vector<double> vec;
  for (const auto& row : Matrix) {
    vec.insert(vec.end(), row.begin(), row.end());
  }
  return vec;
}

/**
 * Converts a 1D vector to a 2D grid data matrix by filling the matrix row by row.
 * This function assumes the total number of elements in the vector is exactly divisible by the specified number of rows.
 *
 * Parameters:
 *   Vec   - A 1D vector containing the grid data, where elements are arranged in row-major order.
 *   NROW  - The number of rows in the resulting matrix.
 *
 * Returns:
 *   A 2D vector (matrix) containing the grid data, arranged by rows.
 */
std::vector<std::vector<double>> GridVec2Mat(const std::vector<double>& Vec,
                                             int NROW){
  // Calculate the number of columns based on the vector size and number of rows
  int NCOL = Vec.size() / NROW;

  // Create the resulting matrix with NROW rows and NCOL columns
  std::vector<std::vector<double>> matrix(NROW, std::vector<double>(NCOL));

  // Fill the matrix with values from the input vector
  for (int i = 0; i < NROW; ++i) {
    for (int j = 0; j < NCOL; ++j) {
      matrix[i][j] = Vec[i * NCOL + j];
    }
  }

  return matrix;
}

/**
 * Computes lagged neighbor values for each cell in a grid using a Moore
 * neighborhood at a specified spatial lag distance, with optional directional
 * filtering.
 *
 * Parameters:
 *   mat    - A 2D vector representing the grid data.
 *   lagNum - The spatial lag distance. A value of 0 returns the current
 *            value of each grid cell.
 *   dir    - Direction selector:
 *              {0}: retain all eight directions.
 *              1: NW, 2: N, 3: NE, 4: W,
 *              5: E, 6: SW, 7: S, 8: SE.
 *            Multiple directions can be specified.
 *
 * Returns:
 *   A 2D vector in which each row corresponds to a grid cell and contains
 *   the lagged values of its selected neighbors. Neighbors outside the
 *   grid are represented by NaN.
 *
 * Notes:
 *   - For lagNum > 0, only cells on the outer Moore ring with Chebyshev
 *     distance equal to lagNum are considered.
 *   - Directional filtering is performed while constructing the neighbor
 *     offsets, avoiding the generation of unnecessary neighbor values.
 *   - Directions are represented internally by an 8-bit mask, with one bit
 *     assigned to each of the eight compass directions.
 *   - The ordering of retained neighbors follows the row-major traversal
 *     of the Moore ring.
 */
std::vector<std::vector<double>> CppLaggedVal4Grid(
    const std::vector<std::vector<double>>& mat,
    int lagNum,
    const std::vector<int>& dir = {0}
) {
  // Validate input
  if (mat.empty() || mat[0].empty() || lagNum < 0) {
    return {};
  }

  const int rows = static_cast<int>(mat.size());
  const int cols = static_cast<int>(mat[0].size());
  const int numCells = rows * cols;

  // lagNum = 0: return the current state of each grid cell
  if (lagNum == 0) {
    std::vector<std::vector<double>> result;
    result.reserve(numCells);

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        result.push_back({mat[i][j]});
      }
    }

    return result;
  }

  // -------------------------------------------------------------------------
  // Construct an 8-bit direction mask.
  //
  // Bit mapping:
  //   bit 0: NW
  //   bit 1: N
  //   bit 2: NE
  //   bit 3: W
  //   bit 4: E
  //   bit 5: SW
  //   bit 6: S
  //   bit 7: SE
  // -------------------------------------------------------------------------

  uint8_t dirMask = 0xFF;  // Default: retain all directions

  if (!(dir.size() == 1 && dir[0] == 0)) {
    dirMask = 0;

    for (int d : dir) {
      if (d >= 1 && d <= 8) {
        dirMask |= static_cast<uint8_t>(1u << (d - 1));
      }
    }
  }

  // -------------------------------------------------------------------------
  // Generate the selected offsets on the outer Moore ring.
  // -------------------------------------------------------------------------

  std::vector<std::pair<int, int>> offsets;
  offsets.reserve(8 * lagNum);

  for (int dx = -lagNum; dx <= lagNum; ++dx) {
    for (int dy = -lagNum; dy <= lagNum; ++dy) {

      // Keep only the outer Moore ring
      if (std::max(std::abs(dx), std::abs(dy)) != lagNum) {
        continue;
      }

      // Directions are assigned by the sign of (row offset, col offset):
      //   1: NW (dx < 0, dy < 0)
      //   2: N  (dx < 0, dy == 0)
      //   3: NE (dx < 0, dy > 0)
      //   4: W  (dx == 0, dy < 0)
      //   5: E  (dx == 0, dy > 0)
      //   6: SW (dx > 0, dy < 0)
      //   7: S  (dx > 0, dy == 0)
      //   8: SE (dx > 0, dy > 0)

      int direction;

      if (dx < 0) {
        direction = (dy < 0) ? 1 : (dy == 0 ? 2 : 3);
      }
      else if (dx == 0) {
        direction = (dy < 0) ? 4 : 5;
      }
      else {
        direction = (dy < 0) ? 6 : (dy == 0 ? 7 : 8);
      }

      // Test the corresponding bit in the direction mask
      if (dirMask & static_cast<uint8_t>(1u << (direction - 1))) {
        offsets.emplace_back(dx, dy);
      }
    }
  }

  const int numNeighbors = static_cast<int>(offsets.size());

  // -------------------------------------------------------------------------
  // Initialize result with NaN.
  // -------------------------------------------------------------------------

  const double NaN = std::numeric_limits<double>::quiet_NaN();

  std::vector<std::vector<double>> result(
    numCells,
    std::vector<double>(numNeighbors, NaN)
  );

  // -------------------------------------------------------------------------
  // Populate lagged neighbor values.
  // -------------------------------------------------------------------------

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {

      const int cellIndex = i * cols + j;

      for (size_t k = 0; k < offsets.size(); ++k) {

        const auto& [dx, dy] = offsets[k];

        const int ni = i + dx;
        const int nj = j + dy;

        if (ni >= 0 && ni < rows &&
            nj >= 0 && nj < cols) {
          result[cellIndex][k] = mat[ni][nj];
        }
      }
    }
  }

  return result;
}

/**
 * Generates grid embeddings by calculating lagged variables for each element in a grid matrix,
 * and stores the results in a matrix where each row represents an element and each column represents
 * a different lagged value or the original element.
 *
 * Parameters:
 *   mat   - A 2D vector representing the grid data.
 *   E     - The number of embedding dimensions (columns in the resulting matrix).
 *   tau   - The spatial lag step for constructing lagged state-space vectors.
 *   style - Embedding style selector:
 *             - style = 0: embedding includes current state as the first dimension.
 *             - style != 0: embedding excludes current state.
 *   dir   - Direction selector (optional):
 *             - If dir = {0}, returns all directional lag values (no filtering).
 *             - If dir ∈ {1,...,8}, keeps only lag values in that direction:
 *                 1: NW, 2: N, 3: NE, 4: W, 5: E, 6: SW, 7: S, 8: SE
 *             - Multiple directions can be specified (e.g., {1,2,3} = NW, N, NE).
 *
 * Returns:
 *   A 2D vector (matrix) where each row contains the averaged lagged variables for
 *   each embedding dimension (column). Columns where all values are NaN are removed.
 *
 * Note:
 *   When tau = 0, lagged variables are calculated for lag steps 0, 1, ..., E-1.
 *   When tau > 0 and style = 0, lagged variables are calculated for lag steps 0, tau, 2*tau, ..., (E-1)*tau.
 *   When tau > 0 and style != 0, lagged variables are calculated for lag steps tau, 2*tau, ..., E*tau.
 */
std::vector<std::vector<double>> GenGridEmbeddings(
    const std::vector<std::vector<double>>& mat,
    int E,
    int tau,
    int style = 1,
    const std::vector<int>& dir = {0}
) {
  int numRows = mat.size();
  int numCols = mat[0].size();
  int total_elements = numRows * numCols;

  // Initialize the result matrix with total_elements rows and E columns and fill with NaN values
  std::vector<std::vector<double>> result(total_elements, std::vector<double>(E, std::numeric_limits<double>::quiet_NaN()));

  if (tau == 0) {
    // Flatten the matrix (mat) into the first column of the result matrix
    int row = 0;
    for (const auto& subset : mat) {
      for (double value : subset) {
        result[row][0] = value;
        ++row;
      }
    }

    // Fill the remaining columns (2 to E) with the averaged lagged variables
    for (int lagNum = 1; lagNum < E; ++lagNum) {
      // Calculate the lagged variables for the current lagNum
      std::vector<std::vector<double>> lagged_vars = CppLaggedVal4Grid(mat, lagNum, dir);

      // Check if all elements in lagged_vars are NaN
      bool allNaN = true;
      for (const auto& subset : lagged_vars) {
        for (double val : subset) {
          if (!std::isnan(val)) {
            allNaN = false;
            break;
          }
        }
        if (!allNaN) break;
      }

      // If all elements are NaN, stop further processing for larger lagNum
      // (since larger lag distances will only produce more NaN results)
      if (allNaN) {
        break;
      }

      // Fill the current column (lagNum) with the averaged lagged variables
      row = 0;
      for (const auto& subset : lagged_vars) {
        double sum = 0.0;
        int count = 0;
        for (double val : subset) {
          if (!std::isnan(val)) {
            sum += val;
            ++count;
          }
        }

        if (count > 0) {
          result[row][lagNum] = sum / count; // Average the valid values
        }
        ++row;
      }
    }
  } else {
    // When tau != 0, calculate lagged variables based on style
    if (style == 0) {
      // style == 0: include current state; lag steps 0, tau, 2*tau, ..., (E-1)*tau
      // Flatten the matrix (mat) into the first column of the result matrix
      int row = 0;
      for (const auto& subset : mat) {
        for (double value : subset) {
          result[row][0] = value;
          ++row;
        }
      }

      // Fill the remaining columns (1 to E-1) with the averaged lagged variables
      for (int i = 1; i < E; ++i) {
        int lagNum = i * tau;

        // Calculate the lagged variables for the current lagNum
        std::vector<std::vector<double>> lagged_vars = CppLaggedVal4Grid(mat, lagNum, dir);

        // Check if all elements in lagged_vars are NaN
        bool allNaN = true;
        for (const auto& subset : lagged_vars) {
          for (double val : subset) {
            if (!std::isnan(val)) {
              allNaN = false;
              break;
            }
          }
          if (!allNaN) break;
        }

        // If all elements are NaN, stop further processing for larger lagNum
        // (since larger lag distances will only produce more NaN results)
        if (allNaN) {
          break;
        }

        // Fill the current column (i) with the averaged lagged variables
        row = 0;
        for (const auto& subset : lagged_vars) {
          double sum = 0.0;
          int count = 0;
          for (double val : subset) {
            if (!std::isnan(val)) {
              sum += val;
              ++count;
            }
          }

          if (count > 0) {
            result[row][i] = sum / count; // Average the valid values
          }
          ++row;
        }
      }
    } else {
      // style != 0: lag steps tau, 2*tau, ..., E*tau (exclude current state)
      for (int i = 1; i <= E; ++i) {
        int lagNum = i * tau;  // Calculate the actual lag step

        // Calculate the lagged variables for the current lagNum
        std::vector<std::vector<double>> lagged_vars = CppLaggedVal4Grid(mat, lagNum, dir);

        // Check if all elements in lagged_vars are NaN
        bool allNaN = true;
        for (const auto& subset : lagged_vars) {
          for (double val : subset) {
            if (!std::isnan(val)) {
              allNaN = false;
              break;
            }
          }
          if (!allNaN) break;
        }

        // If all elements are NaN, stop further processing for larger lagNum
        // (since larger lag distances will only produce more NaN results)
        if (allNaN) {
          break;
        }

        // Fill the current column (i-1) with the averaged lagged variables
        int row = 0;
        for (const auto& subset : lagged_vars) {
          double sum = 0.0;
          int count = 0;
          for (double val : subset) {
            if (!std::isnan(val)) {
              sum += val;
              ++count;
            }
          }

          if (count > 0) {
            result[row][i - 1] = sum / count; // Average the valid values
          }
          ++row;
        }
      }
    }
  }

  // Calculate validColumns (indices of columns that are not entirely NaN)
  std::vector<size_t> validColumns; // To store indices of valid columns

  // Iterate over each column to check if it contains any non-NaN values
  for (size_t col = 0; col < result[0].size(); ++col) {
    bool isAllNaN = true;
    for (size_t row = 0; row < result.size(); ++row) {
      if (!std::isnan(result[row][col])) {
        isAllNaN = false;
        break;
      }
    }
    if (!isAllNaN) {
      validColumns.push_back(col); // Store the index of valid columns
    }
  }

  // If no columns are removed, return the original result
  if (validColumns.size() == result[0].size()) {
    return result;
  } else {
    // // Issue a warning if any columns are removed
    // std::cerr << "Warning: remove all-NA embedding vector columns caused by excessive embedding dimension E selection." << std::endl;

    // Construct the filtered embeddings matrix
    std::vector<std::vector<double>> filteredEmbeddings;
    filteredEmbeddings.reserve(validColumns.size());

    for (size_t row = 0; row < result.size(); ++row) {
      std::vector<double> filteredRow;
      for (size_t col : validColumns) {
        filteredRow.push_back(result[row][col]);
      }
      filteredEmbeddings.push_back(std::move(filteredRow));
    }

    // Return the filtered embeddings matrix
    return filteredEmbeddings;
  }
}

/**
 * Generates multi-level grid embeddings by computing lagged neighbor values
 * for each embedding dimension and returning all lag results instead of averaging them.
 *
 * Each embedding dimension corresponds to a lagged state (τ, 2τ, ... depending on style),
 * and each lagged state is stored as a separate 2D structure representing the grid's lagged values.
 *
 * Parameters:
 *   mat   - A 2D vector representing the grid data.
 *   E     - The number of embedding dimensions (number of lag steps to compute).
 *   tau   - The spatial lag step for constructing lagged state-space vectors.
 *   style - Embedding style selector:
 *             - style = 0: embedding includes the current state as the first dimension.
 *             - style != 0: embedding excludes the current state.
 *   dir   - Direction selector (optional):
 *             - If dir = {0}, returns all directional lag values (no filtering).
 *             - If dir ∈ {1,...,8}, keeps only lag values in that direction:
 *                 1: NW, 2: N, 3: NE, 4: W, 5: E, 6: SW, 7: S, 8: SE
 *             - Multiple directions can be specified (e.g., {1,2,3} = NW, N, NE).
 *
 * Returns:
 *   A 3D vector (std::vector<std::vector<std::vector<double>>>) structured as:
 *     - Outer level (size E): Each element corresponds to an embedding dimension.
 *     - Middle level: Each element corresponds to a grid cell (flattened from 2D).
 *     - Inner level: Contains all lagged neighbor values (may include NaNs).
 *
 * Notes:
 *   - This function differs from GenGridEmbeddings in that it does NOT average lagged values.
 *     Instead, it preserves the full neighbor value sets for each lag step.
 *   - This function also enables direction-based filtering:
 *      - For example, for a 3x3 window (lag = 1), directions map to indices 0–7:
 *         1: NW (0), 2: N (1), 3: NE (2), 4: W (3), 5: E (4), 6: SW (5), 7: S (6), 8: SE (7)
 *      - For larger lags, each direction expands radially; indices are extended accordingly.
 *   - For tau and style:
 *      - When tau = 0, lag steps are sequential (0, 1, 2, ..., E-1).
 *      - When tau > 0 and style = 0, lag steps are 0, τ, 2τ, ..., (E-1)τ.
 *      - When tau > 0 and style != 0, lag steps are τ, 2τ, ..., Eτ.
 */
std::vector<std::vector<std::vector<double>>> GenGridEmbeddingsCom(
    const std::vector<std::vector<double>>& mat,
    int E,
    int tau,
    int style = 1,
    const std::vector<int>& dir = {0}
) {
  std::vector<std::vector<std::vector<double>>> embeddings; // Final 3D result container
  embeddings.reserve(E); // Reserve space for E embedding levels

  if (tau == 0) {
    // tau = 0: lag steps are 0, 1, 2, ..., E-1
    for (int lagNum = 0; lagNum < E; ++lagNum) {
      std::vector<std::vector<double>> lagged_vals = CppLaggedVal4Grid(mat, lagNum, dir);

      // Check if all elements in lagged_vals are NaN
      bool allNaN = true;
      for (const auto& subset : lagged_vals) {
        for (double val : subset) {
          if (!std::isnan(val)) {
            allNaN = false;
            break;
          }
        }
        if (!allNaN) break;
      }

      // If all elements are NaN, stop further processing for larger lagNum
      // (since larger lag distances will only produce more NaN results)
      if (allNaN) {
        break;
      }

      embeddings.push_back(std::move(lagged_vals));
    }
  } else {
    // tau > 0 cases
    if (style == 0) {
      // style == 0: include current state; lag steps: 0, τ, 2τ, ..., (E-1)τ
      for (int i = 0; i < E; ++i) {
        int lagNum = i * tau;
        std::vector<std::vector<double>> lagged_vals = CppLaggedVal4Grid(mat, lagNum, dir);

        // Check if all elements in lagged_vals are NaN
        bool allNaN = true;
        for (const auto& subset : lagged_vals) {
          for (double val : subset) {
            if (!std::isnan(val)) {
              allNaN = false;
              break;
            }
          }
          if (!allNaN) break;
        }

        // If all elements are NaN, stop further processing for larger lagNum
        // (since larger lag distances will only produce more NaN results)
        if (allNaN) {
          break;
        }

        embeddings.push_back(std::move(lagged_vals));
      }
    } else {
      // style != 0: exclude current state; lag steps: τ, 2τ, ..., Eτ
      for (int i = 1; i <= E; ++i) {
        int lagNum = i * tau;
        std::vector<std::vector<double>> lagged_vals = CppLaggedVal4Grid(mat, lagNum, dir);

        // Check if all elements in lagged_vals are NaN
        bool allNaN = true;
        for (const auto& subset : lagged_vals) {
          for (double val : subset) {
            if (!std::isnan(val)) {
              allNaN = false;
              break;
            }
          }
          if (!allNaN) break;
        }

        // If all elements are NaN, stop further processing for larger lagNum
        // (since larger lag distances will only produce more NaN results)
        if (allNaN) {
          break;
        }

        embeddings.push_back(std::move(lagged_vals));
      }
    }
  }

  // // --- Directional filtering section ---
  // if (!(dir.size() == 1 && dir[0] == 0)) {

  //   // Helper function: compute column indices for each of 8 compass directions
  //   // according to the spatial ring pattern around the grid center.
  //   //
  //   // The grid window is (2*lagNum+1) × (2*lagNum+1).
  //   // Only the outer ring (Chebyshev distance == lagNum) is scanned
  //   // in row-major order (top-left to bottom-right).
  //   //
  //   // Directions are assigned by the sign of (row offset, col offset):
  //   //   1: NW (dr < 0, dc < 0)
  //   //   2: N  (dr < 0, dc == 0)
  //   //   3: NE (dr < 0, dc > 0)
  //   //   4: W  (dr == 0, dc < 0)
  //   //   5: E  (dr == 0, dc > 0)
  //   //   6: SW (dr > 0, dc < 0)
  //   //   7: S  (dr > 0, dc == 0)
  //   //   8: SE (dr > 0, dc > 0)
  //   //
  //   // For lagNum > 1, diagonal directions (NW, NE, SW, SE)
  //   // accumulate multiple indices; N/S/E/W remain single.
  //   auto getDirIndicesForLag = [&](int lagNum) {
  //     std::vector<std::vector<int>> dirGroups(8);
  //     int idx = 0; // 0-based column index for C++
  //     for (int dr = -lagNum; dr <= lagNum; ++dr) {
  //       for (int dc = -lagNum; dc <= lagNum; ++dc) {
  //         // only take perimeter cells (Chebyshev distance == lagNum)
  //         if (std::max(std::abs(dr), std::abs(dc)) != lagNum) continue;

  //         // assign this position to its directional group
  //         if (dr < 0 && dc < 0)       dirGroups[0].push_back(idx); // NW
  //         else if (dr < 0 && dc == 0) dirGroups[1].push_back(idx); // N
  //         else if (dr < 0 && dc > 0)  dirGroups[2].push_back(idx); // NE
  //         else if (dr == 0 && dc < 0) dirGroups[3].push_back(idx); // W
  //         else if (dr == 0 && dc > 0) dirGroups[4].push_back(idx); // E
  //         else if (dr > 0 && dc < 0)  dirGroups[5].push_back(idx); // SW
  //         else if (dr > 0 && dc == 0) dirGroups[6].push_back(idx); // S
  //         else if (dr > 0 && dc > 0)  dirGroups[7].push_back(idx); // SE

  //         ++idx;
  //       }
  //     }
  //     return dirGroups;
  //   };

  //   // --- Apply directional filtering for each embedding level ---
  //   for (auto& emb : embeddings) {
  //     if (emb.empty() || emb[0].empty()) continue;

  //     int totalCols = static_cast<int>(emb[0].size());
  //     if (totalCols <= 1) continue; // skip single-column (current-state-only) embeddings

  //     // infer lag number by total column count (8 * lagNum)
  //     int lagNum = totalCols / 8;
  //     if (lagNum < 1) lagNum = 1;

  //     // compute direction index groups for this lag
  //     auto dirGroups = getDirIndicesForLag(lagNum);

  //     // collect all column indices belonging to the requested directions
  //     std::vector<int> selectedCols;
  //     for (int d : dir) {
  //       if (d >= 1 && d <= 8) {
  //         selectedCols.insert(selectedCols.end(),
  //                             dirGroups[d - 1].begin(),
  //                             dirGroups[d - 1].end());
  //       }
  //     }

  //     // ensure unique, sorted indices
  //     std::sort(selectedCols.begin(), selectedCols.end());
  //     selectedCols.erase(std::unique(selectedCols.begin(), selectedCols.end()), selectedCols.end());

  //     // subset each row of this embedding by selected column indices
  //     for (auto& row : emb) {
  //       std::vector<double> filteredRow;
  //       filteredRow.reserve(selectedCols.size());
  //       for (int idx : selectedCols) {
  //         if (idx >= 0 && idx < static_cast<int>(row.size())) {
  //           filteredRow.push_back(row[idx]);
  //         }
  //       }
  //       row = std::move(filteredRow);
  //     }
  //   }
  // }

  // Calculate validSubsets (indices of subsets that are not entirely NaN)
  std::vector<size_t> validSubsets; // To store indices of valid subsets

  // Iterate over each subset to check if it contains any non-NaN values
  for (size_t sub = 0; sub < embeddings.size(); ++sub) {
    bool isAllNaN = true;
    for (size_t row = 0; row < embeddings[sub].size(); ++row) {
      for (size_t col = 0; col < embeddings[sub][row].size(); ++col) {
        if (!std::isnan(embeddings[sub][row][col])) {
          isAllNaN = false;
          break;
        }
      }
      if (!isAllNaN) {
        break;
      }
    }

    if (!isAllNaN) {
      validSubsets.push_back(sub); // Store the index of valid subsets
    }
  }

  // If no subsets are removed, return the original embeddings
  if (validSubsets.size() == embeddings.size()) {
    return embeddings;
  } else {
    // Construct the filtered embeddings
    std::vector<std::vector<std::vector<double>>> filteredEmbeddings;
    filteredEmbeddings.reserve(validSubsets.size());
    for (size_t sub : validSubsets) {
      filteredEmbeddings.push_back(std::move(embeddings[sub]));
    }

    // Return the filtered embeddings
    return filteredEmbeddings;
  }
}

/**
 * @brief Generate k nearest neighbors for all cells in a grid,
 *        choosing only from cells listed in lib and using Queen adjacency.
 *
 * This function identifies, for each cell in a 2D grid (represented by `mat`), its k nearest neighbors
 * among a set of valid cells specified by `lib`. Neighbor relationships are determined based on
 * Queen contiguity (8 directions: N, NE, E, SE, S, SW, W, NW), recursively expanding outward if fewer
 * than k valid neighbors are found in immediate adjacency. The "distance" between cells is calculated
 * as the absolute difference in cell values. When more than k valid candidates are available, the k
 * closest cells (in value) are selected.
 *
 * The output is a vector of size mat.size() × k, where each entry corresponds to the linear index
 * (row-major) of one neighbor per cell, and the output is grouped per row cell.
 *
 * @param mat 2D matrix of values (with possible NaNs).
 * @param lib List of valid cells (flattened indices) to choose neighbors from.
 * @param k   Number of neighbors to find for each cell.
 * @return    A list of neighbors for each cell in mat, in row-major order.
 */
std::vector<std::vector<int>> GenGridNeighbors(
    const std::vector<std::vector<double>>& mat,
    const std::vector<int>& lib,
    size_t k)
{
  size_t nrow = mat.size();
  if (nrow == 0) return {};
  size_t ncol = mat[0].size();
  size_t total_cells = nrow * ncol;

  // Fast lookup for lib
  std::unordered_set<int> lib_set(lib.begin(), lib.end());

  // Queen directions (8 neighbors)
  const std::vector<std::pair<int, int>> directions = {
    {-1, -1}, {-1, 0}, {-1, 1},
    { 0, -1},          { 0, 1},
    { 1, -1}, { 1, 0}, { 1, 1}
  };

  std::vector<std::vector<int>> result(total_cells);  // One entry per cell

  for (size_t idx = 0; idx < total_cells; ++idx) {
    int row = idx / ncol;
    int col = idx % ncol;

    double center_val = mat[row][col];

    std::unordered_set<int> visited;
    visited.insert(static_cast<int>(idx));

    std::set<int> neighbor_candidates;

    std::queue<std::pair<int, int>> frontier;
    frontier.push({row, col});

    while (!frontier.empty() && neighbor_candidates.size() < k) {
      auto [r, c] = frontier.front(); frontier.pop();

      for (const auto& [dr, dc] : directions) {
        int nr = r + dr;
        int nc = c + dc;
        if (nr < 0 || nr >= static_cast<int>(nrow) ||
            nc < 0 || nc >= static_cast<int>(ncol)) continue;

        int flat_nb = nr * ncol + nc;
        if (visited.count(flat_nb)) continue;
        visited.insert(flat_nb);

        // Only consider neighbors that are in lib and not NaN
        if (lib_set.count(flat_nb) && !std::isnan(mat[nr][nc])) {
          neighbor_candidates.insert(flat_nb);
        }

        // Always expand the search
        frontier.push({nr, nc});
      }
    }

    // If more than k, pick those closest in value
    std::vector<std::pair<double, int>> dists;
    for (int nb_idx : neighbor_candidates) {
      int nr = nb_idx / ncol;
      int nc = nb_idx % ncol;
      double dist = std::abs(mat[nr][nc] - center_val);
      dists.emplace_back(dist, nb_idx);
    }

    std::sort(dists.begin(), dists.end());

    std::vector<int> neighbors;
    for (size_t i = 0; i < std::min(k, dists.size()); ++i) {
      neighbors.push_back(dists[i].second);
    }

    result[idx] = std::move(neighbors);
  }

  return result;
}

/**
 * @brief Perform grid-based symbolization on a 2D numeric matrix using local neighborhood statistics.
 *
 * This function computes a symbolic representation (`fs`) for each non-NaN grid cell in `pred` based on its
 * k most similar neighbors selected from the `lib` set. Each symbol reflects local spatial consistency with
 * respect to a global threshold derived from known values.
 *
 * The updated process includes:
 * 1. Compute the global threshold (`s_me`) as the median of all valid (non-NaN) values from `lib` locations only.
 * 2. For each prediction location in `pred`:
 *    - Identify neighbors from `lib` by radially expanding the neighborhood (Queen's case) until at least `k` valid neighbors are found.
 *    - Rank neighbors by absolute difference from the center cell's value.
 *    - Select the top `k` neighbors and assign each a binary symbol (`tau_s`) based on comparison to the global median (`s_me`).
 *    - Compare each `tau_s` to the center cell's own binary symbol (`taus`), and count how many are consistent.
 *    - Store this count as the final symbolic value `fs` for the cell.
 *
 * @param mat A 2D grid (matrix) of values to be symbolized. `NaN` values are treated as missing.
 * @param lib A set of valid reference (library) cell locations as (row, col) pairs. Used for both thresholding and neighbor selection.
 * @param pred A set of prediction cell locations as (row, col) pairs where symbols will be computed.
 * @param k The number of nearest neighbors (by value similarity) to use in symbolization.
 *
 * @return A 1D vector of symbolic values corresponding to the prediction locations, in row-major order.
 *         Entries remain `NaN` if the center cell is invalid or lacks enough valid neighbors.
 *
 * Notes:
 * - Uses Queen’s case (8-directional) neighborhood expansion for neighbor search.
 * - Median is computed only from valid values in the `lib` set (not the whole matrix).
 * - Particularly useful in spatial symbolic analysis, spatial entropy estimation, or causal lattice models.
 */
std::vector<double> GenGridSymbolization(
    const std::vector<std::vector<double>>& mat,
    const std::vector<std::pair<int, int>>& lib,
    const std::vector<std::pair<int, int>>& pred,
    size_t k) {

  // // Validate matrix structure
  // if (mat.empty() || mat[0].empty()) {
  //   throw std::invalid_argument("Invalid matrix dimensions");
  // }
  const size_t rows = mat.size();
  const size_t cols = mat[0].size();

  // Validate and prepare library set
  std::unordered_set<size_t> lib_indices;
  for (const auto& [i,j] : lib) {
    // if (i < 0 || i >= rows || j < 0 || j >= cols) {
    //   throw std::invalid_argument("Library index out of bounds");
    // }
    lib_indices.insert(static_cast<size_t>(i * static_cast<int>(cols) + j));
  }

  // // Validate prediction indices
  // for (const auto& [i,j] : pred) {
  //   if (i < 0 || i >= rows || j < 0 || j >= cols) {
  //     throw std::invalid_argument("Prediction index out of bounds");
  //   }
  // }

  std::vector<double> valid_values;

  // // Compute global median from all non-NaN values
  // for (const auto& row : mat) {
  //   for (double val : row) {
  //     if (!std::isnan(val)) valid_values.push_back(val);
  //   }
  // }

  // Compute the median using only values at lib indices
  for (const auto& [i,j] : lib) {
    double val = mat[i][j];
    if (!std::isnan(val)) valid_values.push_back(val);
  }
  const double s_me = CppMedian(valid_values, true);

  // Prepare result vector
  std::vector<double> result(pred.size(), std::numeric_limits<double>::quiet_NaN());

  // Process each prediction location
  for (size_t idx = 0; idx < pred.size(); ++idx) {
    const auto [i, j] = pred[idx];
    const double center = mat[i][j];
    if (std::isnan(center)) continue;

    // Collect valid neighbors from library
    std::vector<std::pair<double, double>> neighbors;
    const int max_lag = std::max({
      static_cast<int>(i),
      static_cast<int>(rows-1 - i),
      static_cast<int>(j),
      static_cast<int>(cols-1 - j)
    });

    // Radial expansion with library constraints
    for (int lag = 1; lag <= max_lag; ++lag) {
      for (int dx = -lag; dx <= lag; ++dx) {
        for (int dy = -lag; dy <= lag; ++dy) {
          if (std::max(std::abs(dx), std::abs(dy)) != lag) continue;

          const int ni = i + dx;
          const int nj = j + dy;
          if (ni < 0 || ni >= static_cast<int>(rows) ||
              nj < 0 || nj >= static_cast<int>(cols)) continue;

          // Check library membership
          if (!lib_indices.count(static_cast<size_t>(ni * cols + nj))) continue;

          const double val = mat[ni][nj];
          if (!std::isnan(val)) {
            neighbors.emplace_back(std::abs(val - center), val);
          }
        }
      }
      if (neighbors.size() >= k) break;
    }

    // Process if sufficient neighbors found
    if (neighbors.size() >= k) {
      std::sort(neighbors.begin(), neighbors.end(),
                [](const auto& a, const auto& b) {
                  return a.first != b.first ? a.first < b.first : a.second < b.second;
                });

      // Calculate symbolic value
      const double taus = center >= s_me ? 1.0 : 0.0;
      double fs = 0.0;
      for (size_t n = 0; n < k; ++n) {
        const double neighbor_val = neighbors[n].second;
        fs += (neighbor_val >= s_me) == taus ? 1.0 : 0.0;
      }

      result[idx] = fs;
    }
  }

  return result;
}

/**
 * @brief Divide a 2D grid (matrix) into approximately square or shaped blocks.
 *
 * This function partitions a 2D matrix into `b` blocks of roughly equal size.
 * The matrix is represented as a vector of row vectors and assumed to be
 * row-major (i.e., each inner vector represents a row of the matrix).
 *
 * The grid can be divided in different ways based on the `shape` parameter:
 * - 1: Horizontal cuts (blocks are divided row-wise).
 * - 2: Vertical cuts (blocks are divided column-wise).
 * - 3: Radial cuts (blocks are divided in a radial pattern from the center).
 *
 * The grid is divided by first estimating a grid layout of `br` rows and `bc` columns
 * such that br * bc >= b and the blocks are as square as possible for shapes 1 and 2.
 * For radial cuts, the grid is divided into concentric rings and sectors.
 *
 * Each cell in the matrix is assigned a block ID ranging from 0 to b-1, stored in
 * a 1D vector corresponding to the flattened row-major order of the matrix.
 *
 * Any leftover blocks are merged into the last block.
 *
 * @param mat A 2D grid represented as a vector of vectors (row-major).
 * @param b   Number of blocks to divide the grid into.
 * @param shape The shape of the cuts (1: horizontal, 2: vertical, 3: radial).
 * @return A vector of size rows * cols where each element is the block label
 *         assigned to the corresponding cell.
 */
std::vector<int> CppDivideGrid(const std::vector<std::vector<double>>& mat,
                               int b, int shape = 3) {
  int rows = static_cast<int>(mat.size());
  if (rows == 0 || b <= 0) return {};

  int cols = static_cast<int>(mat[0].size());
  int total_cells = rows * cols;

  std::vector<int> labels(total_cells, -1);

  if (shape == 1) {
    // Horizontal cuts: divide the grid row-wise
    int block_h = rows / b;
    int remaining_h = rows % b;

    int row_start = 0;
    for (int i = 0; i < b; ++i) {
      int row_end = row_start + block_h + (i < remaining_h ? 1 : 0);

      for (int r = row_start; r < row_end; ++r) {
        for (int c = 0; c < cols; ++c) {
          int index = r * cols + c;
          labels[index] = i;
        }
      }

      row_start = row_end;
    }
  } else if (shape == 2) {
    // Vertical cuts: divide the grid column-wise
    int block_w = cols / b;
    int remaining_w = cols % b;

    int col_start = 0;
    for (int i = 0; i < b; ++i) {
      int col_end = col_start + block_w + (i < remaining_w ? 1 : 0);

      for (int r = 0; r < rows; ++r) {
        for (int c = col_start; c < col_end; ++c) {
          int index = r * cols + c;
          labels[index] = i;
        }
      }

      col_start = col_end;
    }
  } else if (shape == 3) {
    // Radial cuts: divide the grid into concentric rings and sectors
    double center_x = (cols - 1) / 2.0;
    double center_y = (rows - 1) / 2.0;
    // double max_radius = std::sqrt(center_x * center_x + center_y * center_y);

    int sectors = b;
    double sector_angle = 2 * M_PI / sectors;

    for (int r = 0; r < rows; ++r) {
      for (int c = 0; c < cols; ++c) {
        double dx = c - center_x;
        double dy = r - center_y;
        // double radius = std::sqrt(dx * dx + dy * dy);
        double angle = std::atan2(dy, dx);

        if (angle < 0) angle += 2 * M_PI;

        int sector = static_cast<int>(angle / sector_angle);
        int index = r * cols + c;
        labels[index] = sector % b;
      }
    }
  } else {
    // Default to square blocks
    int br = static_cast<int>(std::sqrt(b));
    int bc = (b + br - 1) / br; // ceiling(b / br)
    while (br * bc < b) ++bc;

    int block_h = rows / br;
    int block_w = cols / bc;
    int remaining_h = rows % br;
    int remaining_w = cols % bc;

    int block_id = 0;

    int row_start = 0;
    for (int i = 0; i < br; ++i) {
      int row_end = row_start + block_h + (i < remaining_h ? 1 : 0);

      int col_start = 0;
      for (int j = 0; j < bc; ++j) {
        int col_end = col_start + block_w + (j < remaining_w ? 1 : 0);

        for (int r = row_start; r < row_end; ++r) {
          for (int c = col_start; c < col_end; ++c) {
            if (r < rows && c < cols) {
              int index = r * cols + c;
              labels[index] = (block_id < b) ? block_id : (b - 1); // fill extras into last block
            }
          }
        }

        ++block_id;
        col_start = col_end;
      }

      row_start = row_end;
    }
  }

  return labels;
}

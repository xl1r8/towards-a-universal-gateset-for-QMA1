#pragma once

#include "eigen.h"
#include "graph.h"
#include "simplex.h"
#include <optional>

using namespace Eigen;

struct CliqueHomology {
  explicit CliqueHomology(Graph g);

  void process();

  const SparseMatrix<double> &diff(int k, bool print = false);
  SparseMatrix<double> &laplacian(int k);

  Graph graph;
  std::vector<std::unordered_map<Simplex, int>> simplex_to_id;

  void print_simplex(std::ostream &stream, const Simplex &s) const;

  int euler_characteristic() const;

private:
  std::array<std::optional<SparseMatrix<double>>, SIMPLEX_MAX_SIZE + 2> diff_cache;
  std::array<std::optional<SparseMatrix<double>>, SIMPLEX_MAX_SIZE + 2> laplacian_cache;
};
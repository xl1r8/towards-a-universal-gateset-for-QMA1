#pragma once

#include "graph.h"
#include "homology.h"
#include "simplex.h"
#include <memory>

constexpr std::string v0 = "v0";
constexpr std::string vertex_names[] = {"xx", "a2", "a3", "a4", "b2", "b3", "b4"};

struct VertexNames {
  VertexNames(int n, const std::string &suffix = "");

  [[nodiscard]] StrVecVec bit_simplices(bool bit, int i) const;

  int n;
  StrVec x, a2, a3, a4, b2, b3, b4;
  std::vector<std::array<std::string, 7>> qubit_vertices;
};

struct QubitGraph : public VertexNames {
  explicit QubitGraph(int n);

  void make_homology();

  [[nodiscard]] VectorXd simplices_to_vec(const StrVecVec &simplices) const;

  [[nodiscard]] VectorXd embed_bits(const BoolVec &bits) const;

  [[nodiscard]] StrVecVec cliques_for_states(const BoolVecVec &states) const;

  [[nodiscard]] VectorXd embed_qubits(const std::vector<std::array<double, 2>> &bits) const;

  [[nodiscard]] int gadget_dim() const;

  void fill_cycle(const Graph &cycle, const StrStrMap &f = {}, const std::string &v0 = "v0");

  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_flip2(bool b0, bool b1, bool minus = true) const;

  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_00p10p11() const;

  [[nodiscard]] Graph make_part2(const StrVecVec &simpl0, const StrVecVec &simpl1, const StrVec &dx,
                                 const StrVec &loop) const;

  //     B = [qket(0) + qket(1), qket(0) - qket(1)]
  //     tq = lambda t, b: tp(tket(t), qket(b), I2)
  //     tqH = lambda t, b: tp(tket(t), B[b], I2)
  //     H_trans = (
  //             kb(tq(0, 0) - tqH(1, 0)) +
  //             kb(tq(0, 1) - tqH(1, 1))
  //     )
  //     tq = lambda t, b: tp(tket(t), I2, qket(b))
  //     tqH = lambda t, b: tp(tket(t), I2, B[b])
  //     H_trans += (
  //             kb(2*tq(1, 0) - tqH(2, 0)) +
  //             kb(2*tq(1, 1) - tqH(2, 1))
  //     )
  //     H = H_trans

  // |00> - |10> - |11>
  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_Htrans_00m10m11() const;

  // |00> - |11> - |11>
  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_00m11m11() const;

  // |01> - |10> + |11>
  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_Htrans_01m10p11() const;

  // |00> + 3|11> [Section 8.3.7, KK23]
  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_00p11p11p11() const;

  [[nodiscard]] std::pair<Graph, StrStrMap> cycle_0m1() const;

  // 2*|10> - |00> - |01>
  // 2*|11> - |00> + |01>

  Graph g;
  std::unique_ptr<CliqueHomology> h;
};
